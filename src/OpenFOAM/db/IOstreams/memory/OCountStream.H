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
    Foam::OSCountStream

Description
    An output stream for calculating byte counts.

\*---------------------------------------------------------------------------*/

#ifndef OScountStream_H
#define OScountStream_H

#include "OSstream.H"
#include <iostream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class countstreambuf Declaration
\*---------------------------------------------------------------------------*/

//- A streambuf class for determining byte counts
class countstreambuf
:
    public std::streambuf
{
    // Private data

        //- The number of bytes
        std::streamsize n_;

protected:
    // Protected members

        //- Put sequence of characters
        virtual std::streamsize xsputn(const char* s, std::streamsize num)
        {
            n_ += num;
            return num;
        }


        //- Set position pointer to absolute position
        //  For the counter, any positioning is ignored and it always acts like
        //  seekpos(0), which resets the count.
        virtual std::streampos seekpos
        (
            std::streampos sp,
            std::ios_base::openmode which = std::ios_base::in|std::ios_base::out
        )
        {
            n_ = 0;
            return 0;
        }


public:

    // Constructors

        //- Construct null
        countstreambuf()
        :
            n_(0)
        {}


    // Access

        //- Get number of bytes counted
        std::streamsize size() const
        {
            return n_;
        }
};


/*---------------------------------------------------------------------------*\
                        Class ocountstream Declaration
\*---------------------------------------------------------------------------*/

//- Trivial output stream for calculating byte counts
//  Since all output values are discarded, it can also be used as a /dev/null
//  output buffer as well
class ocountstream
:
    virtual public std::ios,
    public std::ostream
{
    // Private data

    countstreambuf buf_;

public:

    // Constructors

        //- Construct null
        ocountstream()
        :
            std::ostream(&buf_)
        {}


    // Member Functions

      // Access

        //- This hides both signatures of std::basic_ios::rdbuf()
        countstreambuf* rdbuf()
        {
            return &buf_;
        }


        //- Get number of bytes counted
        std::streamsize size() const
        {
            return buf_.size();
        }


        //- Rewind the stream, reset the count
        void rewind()
        {
            buf_.pubseekpos(0);
            clear();    // for safety, clear any old errors
        }
};


/*---------------------------------------------------------------------------*\
                        Class OCountStream Declaration
\*---------------------------------------------------------------------------*/

//- An output stream for calculating byte counts
class OCountStream
:
    private ocountstream,
    public OSstream
{

    // Private Member Functions

        //- Disallow default bitwise copy construct
        OCountStream(const OCountStream&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const OCountStream&) = delete;

public:

    // Constructors

        //- Construct and set stream status
        OCountStream
        (
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            ocountstream(),
            OSstream
            (
                static_cast<ocountstream&>(*this),
                "output",
                format,
                version
        )
        {}


    //- Destructor
    ~OCountStream()
    {}


    // Member functions

      // Access

        //- Return the number of bytes counted
        using ocountstream::size;


      // Edit

        //- Rewind the stream, reset the count, clearing any old errors
        void rewind()
        {
            ocountstream::rewind();
            setGood();  // resynchronize with internal state
        }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
