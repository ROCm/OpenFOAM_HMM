/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2018 OpenCFD Ltd.
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
            buffer.resize(os.size());       // synchronize sizes
        }

        something.write(buffer, buffer.size());
    \endcode

    Although the UOListStream is quite lightweight, there may be cases
    where it is preferable to reuse the stream as well.
    \code
        DynamicList<char> buffer(4096);     // allocate some large buffer

        UOListStream os(buffer);
        os << "content1" << " and more content";
        buffer.resize(os.size());           // synchronize sizes

        something.write(buffer, buffer.size());

        os.rewind();
        os << "content2";
        buffer.resize(os.size());      // synchronize sizes

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

namespace Detail
{

/*---------------------------------------------------------------------------*\
                Class Detail::UOListStreamAllocator Declaration
\*---------------------------------------------------------------------------*/

//- An stream/stream-buffer allocator for external buffers
class UOListStreamAllocator
{
protected:

    // Protected data

        typedef std::ostream stream_type;

        //- The stream buffer
        memorybuf::out buf_;

        //- The stream
        stream_type stream_;


    // Constructors

        //- Construct for character array and number of bytes
        UOListStreamAllocator(char *buffer, size_t nbytes)
        :
            buf_(buffer, nbytes),
            stream_(&buf_)
        {}

        void printBufInfo(Ostream& os) const
        {
            buf_.printBufInfo(os);
        }

public:

    // Public Member Functions

        //- Const UList access to the characters written (shallow copy).
        inline const UList<char> list() const
        {
            return buf_.list();
        }

        //- Non-const UList access to the characters written (shallow copy).
        inline UList<char> list()
        {
            return buf_.list();
        }

        //- The current list output capacity
        inline label capacity() const
        {
            return buf_.capacity();
        }

        //- The current output position in the buffer,
        //- which is also the addressed list size
        inline label size() const
        {
            return buf_.tellp();
        }

        //- Move to buffer start, clear errors
        void rewind()
        {
            buf_.pubseekpos(0, std::ios_base::out);
            stream_.clear(); // for safety, clear any old errors
        }
};

} // End namespace Detail


/*---------------------------------------------------------------------------*\
                         Class UOListStream Declaration
\*---------------------------------------------------------------------------*/

//- An OSstream attached to an unallocated external buffer
class UOListStream
:
    public Detail::UOListStreamAllocator,
    public OSstream
{
    typedef Detail::UOListStreamAllocator allocator_type;

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
        template<unsigned N>
        UOListStream
        (
            FixedList<char, N>& buffer,
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            UOListStream(buffer.data(), N, format, version)
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


    // Member Functions

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
