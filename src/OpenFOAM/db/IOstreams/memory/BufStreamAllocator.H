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
    Foam::BufStreamAllocator

Description
    Helper for memory buffer streams such as UIBufStream, UOBufStream

\*---------------------------------------------------------------------------*/

#ifndef BufStreamAllocator_H
#define BufStreamAllocator_H

#include <type_traits>
#include <sstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class memorybuf Declaration
\*---------------------------------------------------------------------------*/

//- A streambuf class for using externally allocated memory for its buffer
class memorybuf
:
    public std::streambuf
{
protected:

    // Protected members

        //- Set position pointer to relative position
        virtual std::streampos seekoff
        (
            std::streamoff off,
            std::ios_base::seekdir way,
            std::ios_base::openmode which = std::ios_base::in|std::ios_base::out
        )
        {
            const bool testin  = which & std::ios_base::in;
            const bool testout = which & std::ios_base::out;

            if (way == std::ios_base::beg)
            {
                if (testin)
                {
                    setg(eback(), eback(), egptr());
                    gbump(off);
                }
                if (testout)
                {
                    setp(pbase(), epptr());
                    pbump(off);
                }

                return off;
            }

            if (way == std::ios_base::cur)
            {
                if (testin)
                {
                    gbump(off);
                }
                if (testout)
                {
                    pbump(off);
                }
            }
            else if (way == std::ios_base::end)
            {
                if (testin)
                {
                    gbump(off);
                }
                if (testout)
                {
                    pbump(off);
                }
            }

            if (testin)
            {
                return gptr() - eback();
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
            std::ios_base::openmode which = std::ios_base::in|std::ios_base::out
        )
        {
            return seekoff(pos, std::ios_base::beg, which);
        }


        //- Get sequence of characters
        virtual std::streamsize xsgetn(char* s, std::streamsize n)
        {
            std::streamsize count = 0;

            // some optimization could be possible here
            while (count < n && gptr() < egptr())
            {
                *(s + count++) = *(gptr());
                gbump(1);
            }

            return count;
        }


        //- Put sequence of characters
        virtual std::streamsize xsputn(const char* s, std::streamsize n)
        {
            std::streamsize count = 0;

            // some optimization could be possible here
            while (count < n && pptr() < epptr())
            {
                *(pptr()) = *(s + count++);
                pbump(1);
            }

            return count;
        }


public:

    // Constructors

        //- Construct for specified buffer
        memorybuf(char* buffer, std::streamsize num)
        {
            setg(buffer, buffer, buffer + num);
            setp(buffer, buffer + num);
        }

};


/*---------------------------------------------------------------------------*\
                     Class BufStreamAllocator Declaration
\*---------------------------------------------------------------------------*/

//- An stream/stream-buffer allocator for external buffers
template<class StreamType, bool Manage=false>
class BufStreamAllocator
{
    // Private data

        //- Storage
        char *storage_;

        //- The number of bytes in the storage
        std::streamsize len_;

        //- Reference to the underlying buffer
        memorybuf buf_;

protected:

    // Protected data

        typedef StreamType stream_type;

        //- The stream pointer
        stream_type stream_;


    // Constructors

        //- Construct with buffer and number of bytes
        BufStreamAllocator(char *buffer, size_t nbytes)
        :
            storage_(buffer),
            len_(nbytes),
            buf_(storage_, len_),
            stream_(&buf_)
        {}


    //- Destructor
    ~BufStreamAllocator()
    {
        // Possible cleanup of storage
        if (Manage && storage_)
        {
            delete storage_;
            storage_ = nullptr;
        }
    }


    // Protected Member Functions


        //- Position of the get buffer
        std::streampos tellg() const
        {
            return const_cast<stream_type&>(stream_).tellg();
        }

        //- Position of the put buffer
        std::streampos tellp() const
        {
            return const_cast<stream_type&>(stream_).tellp();
        }


public:

    // Public Member Functions

        //- Move to buffer start, clear errors
        void rewind()
        {
            stream_.rdbuf()->pubseekpos(0);
            stream_.clear(); // for safety, clear any old errors
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
