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
    Foam::memorybuf

Description
    A std::streambuf used for memory buffer streams such as
    UIListStream, UOListStream, etc.

\*---------------------------------------------------------------------------*/

#ifndef memoryStreamBuffer_H
#define memoryStreamBuffer_H

#include "UList.H"
#include <type_traits>
#include <sstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class memorybuf Declaration
\*---------------------------------------------------------------------------*/

//- A streambuf for memory
class memorybuf
:
    public std::streambuf
{
protected:

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
        }
        else if (way == std::ios_base::cur)
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
                setg(eback(), eback(), egptr());
                gbump(egptr() - eback() - off);
            }
            if (testout)
            {
                setp(pbase(), epptr());
                pbump(epptr() - pbase() - off);
            }
        }

        if (testin)
        {
            return (gptr() - eback()); // tellg()
        }
        if (testout)
        {
            return (pptr() - pbase()); // tellp()
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


public:

    // Forward declarations
    class in;
    class out;
};


/*---------------------------------------------------------------------------*\
                        Class memorybuf::in Declaration
\*---------------------------------------------------------------------------*/

//- An input streambuf for memory access
class memorybuf::in
:
    public memorybuf
{
protected:

    //- Get sequence of characters
    virtual std::streamsize xsgetn(char* s, std::streamsize n)
    {
        std::streamsize count = 0;

        while (count < n && gptr() < egptr())
        {
            *(s + count++) = *(gptr());
            gbump(1);
        }

        return count;
    }


    //- The buffer get position
    inline std::streamsize tellg() const
    {
        return (gptr() - eback());
    }

    //- Sync get buffer pointers to agree with list dimensions
    inline void set_gbuffer(UList<char>& list)
    {
        setg(list.begin(), list.begin(), list.end());
    }
};


/*---------------------------------------------------------------------------*\
                       Class memorybuf::out Declaration
\*---------------------------------------------------------------------------*/

//- An output streambuf for memory access
class memorybuf::out
:
    public memorybuf
{
protected:

    //- Put sequence of characters
    virtual std::streamsize xsputn(const char* s, std::streamsize n)
    {
        std::streamsize count = 0;
        while (count < n && pptr() < epptr())
        {
            *(pptr()) = *(s + count++);
            pbump(1);
        }

        return count;
    }

    //- The buffer put position
    inline std::streamsize tellp() const
    {
        return (pptr() - pbase());
    }

    //- Sync put buffer pointers to agree with list dimensions
    inline void set_pbuffer(UList<char>& list)
    {
        setp(list.begin(), list.end());
    }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
