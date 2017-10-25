/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::OListStream

Description
    An output stream that writes to a DynamicList.

\*---------------------------------------------------------------------------*/

#ifndef OListStream_H
#define OListStream_H

#include <sstream>
#include "DynamicList.H"
#include "OSstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration
class OListStreamAllocator;

/*---------------------------------------------------------------------------*\
                    Class OListStreamAllocator Declaration
\*---------------------------------------------------------------------------*/

//- An stream/stream-buffer allocator for external buffers
class OListStreamAllocator
{
    //- A streambuf adapter to output to a DynamicList
    class olistbuf
    :
        public std::streambuf
    {
        friend OListStreamAllocator;

        //- Underlying list storage
        DynamicList<char,512> storage_;


        //- Adjust buffer pointers to agree with list sizes
        inline void syncBufferPointers()
        {
            setp(storage_.data(), storage_.data() + storage_.capacity());
            pbump(storage_.size());
        }

        //- Adjust addressed list size to agree with buffer pointers
        inline void syncListSize()
        {
            storage_.setSize(pptr() - pbase());
        }


    protected:

        // Protected members

            //- Set position pointer to relative position
            virtual std::streampos seekoff
            (
                std::streamoff off,
                std::ios_base::seekdir way,
                std::ios_base::openmode which = std::ios_base::out
            )
            {
                const bool testout = which & std::ios_base::out;

                if (way == std::ios_base::beg)
                {
                    if (testout)
                    {
                        setp(pbase(), epptr());
                        pbump(off);
                    }

                    return off;
                }

                if (way == std::ios_base::cur)
                {
                    if (testout)
                    {
                        pbump(off);
                    }
                }
                else if (way == std::ios_base::end)
                {
                    if (testout)
                    {
                        pbump(off);
                    }
                }

                if (testout)
                {
                    return pptr() - pbase();
                }

                return -1;
            }


            //- Set position pointer to absolute position
            virtual std::streampos seekpos
            (
                std::streampos pos,
                std::ios_base::openmode which = std::ios_base::out
            )
            {
                if (which & std::ios_base::out)
                {
                    setp(pbase(), epptr());
                    pbump(pos);

                    return pptr() - pbase();
                }

                return -1;
            }


            //- Put sequence of characters
            virtual std::streamsize xsputn(const char* s, std::streamsize n)
            {
                const std::streamsize newlen = n + storage_.size();

                if (newlen > storage_.capacity())
                {
                    // Either use reserve(), or setCapacity() directly
                    // with finer control over growth
                    storage_.reserve(newlen);
                    syncBufferPointers();
                }

                std::streamsize count = 0;
                while (count < n && pptr() < epptr())
                {
                    *(pptr()) = *(s + count++);
                    pbump(1);
                }

                // Synchronize list size with output
                syncListSize();

                return count;
            }


    public:

        // Constructors

            //- Construct with an empty list
            olistbuf()
            :
                storage_(1024)
            {
                syncBufferPointers();
            }


            //- Construct with a specified number of reserved bytes
            olistbuf(size_t nbytes)
            :
                storage_()
            {
                storage_.reserve(std::max(label(nbytes),1024));
                syncBufferPointers();
            }


            //- Move construct from an existing List
            olistbuf(List<char>&& buffer)
            :
                storage_(std::move(buffer))
            {
                syncBufferPointers();
            }


            //- Move construct from an existing DynamicList
            template<int AnySize>
            olistbuf(DynamicList<char,AnySize>&& buffer)
            :
                storage_(std::move(buffer))
            {
                syncBufferPointers();
            }
    };


    // Private data

        //- Reference to the underlying buffer
        olistbuf buf_;


protected:

    // Protected data

        typedef std::ostream stream_type;

        //- The output stream
        stream_type stream_;


    // Constructors

        //- Construct with an empty list
        OListStreamAllocator()
        :
            buf_(),
            stream_(&buf_)
        {}


        //- Construct with a specified number of reserved bytes
        OListStreamAllocator(size_t nbytes)
        :
            buf_(nbytes),
            stream_(&buf_)
        {}


        //- Move construct from an existing List
        OListStreamAllocator(List<char>&& buffer)
        :
            buf_(std::move(buffer)),
            stream_(&buf_)
        {}


        //- Move construct from an existing DynamicList
        template<int SizeMin>
        OListStreamAllocator(DynamicList<char,SizeMin>&& buffer)
        :
            buf_(std::move(buffer)),
            stream_(&buf_)
        {}


    //- Destructor
    ~OListStreamAllocator()
    {}


public:

    // Member Functions

        //- Content as a list of characters
        UList<char> list() const
        {
            return UList<char>
            (
                const_cast<char*>(buf_.storage_.cdata()),
                buf_.storage_.size()
            );
        }


        //- Content as a list of characters
        UList<char>& list()
        {
            return static_cast<UList<char>&>(buf_.storage_);
        }

        //- Return the current list capacity
        inline label capacity() const
        {
            return buf_.storage_.capacity();
        }


        //- Reserve allocation space for at least this size.
        inline void reserve(const label nElem)
        {
            buf_.storage_.reserve(nElem);
            buf_.syncBufferPointers();
        }


        //- Return the current output position in the buffer
        //  The same as the DynamicList::size()
        std::streampos size() const
        {
            return const_cast<stream_type&>(stream_).tellp();
        }


        //- Move to buffer start, clear errors
        void rewind()
        {
            buf_.storage_.clear();
            buf_.syncBufferPointers();
            stream_.clear(); // for safety, clear any old errors
        }
};


/*---------------------------------------------------------------------------*\
                         Class OListStream Declaration
\*----------------------------------------------d-----------------------------*/

//- An OSstream attached a DynamicList
class OListStream
:
    public OListStreamAllocator,
    public OSstream
{
    typedef OListStreamAllocator allocator_type;

public:

    // Constructors

        //- Construct with an empty list
        OListStream
        (
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            allocator_type(),
            OSstream(stream_, "output", format,version)
        {}


        //- Construct with a specified number of reserved bytes
        OListStream
        (
            size_t nbytes,
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            allocator_type(nbytes),
            OSstream(stream_, "output", format, version)
        {}


        //- Move construct from an existing List
        OListStream
        (
            List<char>&& buffer,
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            allocator_type(std::move(buffer)),
            OSstream(stream_, "output", format, version)
        {}


        //- Move construct from an existing DynamicList
        template<int SizeMin>
        OListStream
        (
            DynamicList<char,SizeMin>&& buffer,
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            allocator_type(std::move(buffer)),
            OSstream(stream_, "output", format, version)
        {}


    //- Destructor
    ~OListStream()
    {}


    // Member functions

        //- Rewind the stream, clearing any old errors
        using allocator_type::rewind;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
