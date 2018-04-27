/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2018 OpenCFD Ltd.
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
    An output stream that writes to a List and manages the List storage.
    Similar to OStringStream but with a List for its storage instead of
    as string to allow reuse of List contents without copying.

    The default list size is 512-bytes with a 256-byte block increment.
    These values can be changed after construction using the reserve() and
    the setBlockSize() methods.

See Also
    Foam::IListStream
    Foam::UOListStream
    Foam::UIListStream

\*---------------------------------------------------------------------------*/

#ifndef OListStream_H
#define OListStream_H

#include "DynamicList.H"
#include "OSstream.H"
#include "stdFoam.H"
#include "memoryStreamBuffer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace Detail
{

/*---------------------------------------------------------------------------*\
                Class Detail::OListStreamAllocator Declaration
\*---------------------------------------------------------------------------*/

//- An stream/stream-buffer output allocator with DynamicList-like storage
class OListStreamAllocator
{
    //- A streambuf adapter with resizing similar to DynamicList
    class dynbuf
    :
        public memorybuf::out
    {
        friend class OListStreamAllocator;

        //- Helper for block size - small list minimum of 64 bytes.
        constexpr static int min_size(int n)
        {
            return stdFoam::max(64, n);
        }

        //- Block size when resizing the list
        int block_ = 256;

        //- Underlying list storage.
        //  Internally manage like a DynamicList, with its capacity known
        //  from the list size and the addressable size known through the
        //  stream pointers.
        List<char> storage_;

    protected:

        //- Increment capacity directly and adjust buffer pointers to
        //- correspond with the storage size.
        inline void minCapacity
        (
            const std::streamsize n,
            const std::streamsize cur = 0
        )
        {
            const auto newEnd = n + cur;
            if (newEnd > storage_.size())
            {
                auto newCapacity =
                (
                    (storage_.size() + block_)
                  - (storage_.size() % block_)
                );

                while (newCapacity < newEnd)
                {
                    newCapacity += block_;
                }

                // Info<<"request:" << newEnd
                //     << " cur cap:" << storage_.size()
                //     << " new cap:" << newCapacity
                //     << " pos:" << cur
                //     << " incr:" << incr << endl;

                storage_.resize(newCapacity);
                sync_pbuffer();
                pbump(cur);
            }
        }

        //- Define new increment
        inline void setBlockSize(const int i)
        {
            const auto prev = block_;
            block_ = min_size(i);

            if (block_ > prev)
            {
                minCapacity(0, tellp());
            }
        }

        //- Handle overflow
        virtual int overflow(int c = EOF)
        {
            if (c != EOF)
            {
                // Need another output block
                minCapacity(block_, tellp());

                *(pptr()) = c;
                pbump(1);
            }
            return c;
        }


        //- Put sequence of characters
        virtual std::streamsize xsputn(const char* s, std::streamsize n)
        {
            // Enough space so that appends work without problem
            minCapacity(n, tellp());

            std::streamsize count = 0;
            while (count < n && pptr() < epptr())
            {
                *(pptr()) = *(s + count++);
                pbump(1);
            }

            return count;
        }

        //- Initialize put buffer
        void init_pbuffer(const std::streamsize n)
        {
            sync_pbuffer();
            minCapacity(n);
        }


    public:

    // Constructors

        //- Construct with an empty list, or specified number of reserved bytes
        dynbuf(size_t nbytes = 512)
        :
            storage_()
        {
            init_pbuffer(min_size(nbytes));
        }

        //- Move construct from List
        dynbuf(List<char>&& buffer)
        :
            storage_(std::move(buffer))
        {
            init_pbuffer(block_);
        }

        //- Move construct from DynamicList
        template<int AnySize>
        dynbuf(DynamicList<char,AnySize>&& buffer)
        :
            storage_(std::move(buffer))
        {
            init_pbuffer(block_);
        }


    // Member Functions

        //- Return the current list output capacity
        inline label capacity() const
        {
            return storage_.size();
        }

        //- Sync put buffer pointers to agree with list dimensions
        inline void sync_pbuffer()
        {
            resetp(storage_.data(), storage_.size());
        }

        //- Clear storage
        inline void clearStorage()
        {
            storage_.clear();
            sync_pbuffer();
        }

        //- Shrink to addressed storage
        inline void shrink()
        {
            const auto cur = tellp();
            storage_.resize(cur);
            sync_pbuffer();
            pbump(cur);
        }

        //- Transfer list contents to other list
        inline void swap(List<char>& other)
        {
            storage_.resize(tellp());   // Addressed area
            storage_.swap(other);
            sync_pbuffer();
        }
    };

protected:

    // Protected data

        typedef std::ostream stream_type;

        //- The stream buffer
        dynbuf buf_;

        //- The stream
        stream_type stream_;


    // Constructors

        //- Construct with an empty list or specified number of reserved bytes
        OListStreamAllocator(size_t nbytes = 512)
        :
            buf_(nbytes),
            stream_(&buf_)
        {}

        //- Move construct from List
        OListStreamAllocator(List<char>&& buffer)
        :
            buf_(std::move(buffer)),
            stream_(&buf_)
        {}

        //- Move construct from DynamicList
        template<int SizeMin>
        OListStreamAllocator(DynamicList<char,SizeMin>&& buffer)
        :
            buf_(std::move(buffer)),
            stream_(&buf_)
        {}


    // Protected Member Functions

        void printBufInfo(Ostream& os) const
        {
            os  << "put=" << buf_.tellp()
                << "/" << buf_.capacity()
                << " block=" << buf_.block_;
        }

public:

    // Member Functions

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

        //- Reserve output space for at least this amount.
        inline void reserve(const std::streamsize n)
        {
            // Also maintain current position when resizing
            const auto cur = buf_.tellp();
            if (n > cur)
            {
                buf_.minCapacity(n - cur, cur);
            }
        }

        //- Adjust block size for output
        inline void setBlockSize(int n)
        {
            return buf_.setBlockSize(n);
        }

        //- Transfer list contents to other list
        inline void swap(List<char>& other)
        {
            buf_.swap(other);
        }

        //- Shrink to addressed space, should not affect stream.
        inline void shrink()
        {
            buf_.shrink();
        }

        //- Clear storage
        void clearStorage()
        {
            buf_.clearStorage();
            stream_.clear(); // for safety, clear any old errors
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
                         Class OListStream Declaration
\*---------------------------------------------------------------------------*/

//- An OSstream attached to a List
class OListStream
:
    public Detail::OListStreamAllocator,
    public OSstream
{
    typedef Detail::OListStreamAllocator allocator_type;

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
            OSstream(stream_, "output", format, version)
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
