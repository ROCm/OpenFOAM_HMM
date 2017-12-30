/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::UOListStream

Description
    Similar to OStringStream but using an externally managed buffer for
    its output.

    This allows the output buffer to be reused and can make it easier when
    writing out data.  It is the user's responsibility to ensure proper
    synchronization in the sizes. Provided that the external buffer is large
    enough that overflow does not occur, the following usage pattern
    works.

    \code
        DynamicList<char> buffer(4096);     // allocate some large buffer

        {
            UOListStream os(buffer);
            os << "content1" << " and more content";
            buffer.setSize(os.size());      // synchronize sizes
        }

        something.write(buffer, buffer.size());
    \endcode

    Although the UOListStream is quite lightweight, there may be cases
    where it is preferable to reuse the stream as well.
    \code
        DynamicList<char> buffer(4096);     // allocate some large buffer

        UOListStream os(buffer);
        os << "content1" << " and more content";
        buffer.setSize(os.size());          // synchronize sizes

        something.write(buffer, buffer.size());

        os.rewind();
        os << "content2";
        buffer.setSize(os.size());      // synchronize sizes

        something.write(buffer, buffer.size());

        // or simply using the output size directly (without sync)
        os.rewind();
        os << "content3";

        something.write(buffer, os.size());
    \endcode

See Also
    Foam::IListStream
    Foam::OListStream
    Foam::UIListStream

\*---------------------------------------------------------------------------*/

#ifndef UOListStream_H
#define UOListStream_H

#include "DynamicList.H"
#include "FixedList.H"
#include "OSstream.H"
#include "memoryStreamBuffer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class UOListAllocator Declaration
\*---------------------------------------------------------------------------*/

//- An stream/stream-buffer allocator for external buffers
class UOListAllocator
{
    //- A streambuf adapter for output to UList or equivalent
    class ulistbuf
    :
        public memorybuf::out
    {
        friend UOListAllocator;

        //- Underlying list storage
        UList<char> list_;

    public:

      // Constructors

        //- Construct for specified buffer
        ulistbuf(char* buffer, std::streamsize num)
        :
            list_(buffer, num)
        {
            set_pbuffer(list_);
        }


      // Member Functions

        //- Return the current list output capacity
        inline label capacity() const
        {
            return list_.size();
        }

        //- Sync put buffer pointers to agree with list dimensions
        inline void sync_pbuffer_to_list()
        {
            set_pbuffer(list_);
        }
    };


protected:

    // Protected data

        typedef std::ostream stream_type;

        //- The stream buffer
        ulistbuf buf_;

        //- The stream
        stream_type stream_;


    // Constructors

        //- Construct with buffer and number of bytes
        UOListAllocator(char *buffer, size_t nbytes)
        :
            buf_(buffer, nbytes),
            stream_(&buf_)
        {}


    // Protected Member Functions

        void printBufInfo(Ostream& os) const
        {
            os  << "pos=" << buf_.tellp()
                << " capacity=" << buf_.capacity();
        }

public:

    // Public Member Functions

        //- The current list output capacity
        inline label capacity() const
        {
            return buf_.capacity();
        }


        //- The current output position in the buffer,
        //  which is also the addressed list size
        inline label size() const
        {
            return buf_.tellp();
        }


        //- Const access to written contents as a list of characters
        inline const UList<char> list() const
        {
            return UList<char>
            (
                const_cast<char*>(buf_.list_.cdata()),
                buf_.tellp()
            );
        }


        //- Non-const access to written contents as a list of characters
        inline UList<char> list()
        {
            return UList<char>
            (
                const_cast<char*>(buf_.list_.cdata()),
                buf_.tellp()
            );
        }


        //- Move to buffer start, clear errors
        void rewind()
        {
            buf_.pubseekpos(0, std::ios_base::out);
            stream_.clear(); // for safety, clear any old errors
        }
};


/*---------------------------------------------------------------------------*\
                         Class UOListStream Declaration
\*---------------------------------------------------------------------------*/

//- An OSstream attached to an unallocated external buffer
class UOListStream
:
    public UOListAllocator,
    public OSstream
{
    typedef UOListAllocator allocator_type;

public:

    // Constructors

        //- Construct using specified buffer and number of bytes
        UOListStream
        (
            char* buffer,
            size_t nbytes,
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            allocator_type(buffer, nbytes),
            OSstream(stream_, "output", format, version)
        {}


        //- Construct using data area from a List and number of bytes
        UOListStream
        (
            UList<char>& buffer,
            size_t size,
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            UOListStream(buffer.data(), size, format, version)
        {}

        //- Construct using data area from a FixedList
        template<unsigned FixedSize>
        UOListStream
        (
            FixedList<char,FixedSize>& buffer,
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            UOListStream(buffer.data(), FixedSize, format, version)
        {}

        //- Construct using data area from a List and its inherent storage size
        UOListStream
        (
            UList<char>& buffer,
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            UOListStream(buffer.data(), buffer.size(), format, version)
        {}


        //- Construct using data area from a DynamicList and its capacity
        template<int SizeMin>
        UOListStream
        (
            DynamicList<char,SizeMin>& buffer,
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            UOListStream(buffer.data(), buffer.capacity(), format, version)
        {}


    // Member functions

        //- Rewind the stream, clearing any old errors
        virtual void rewind()
        {
            allocator_type::rewind();
            setGood();  // resynchronize with internal state
        }

        //- Print description to Ostream
        virtual void print(Ostream& os) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
