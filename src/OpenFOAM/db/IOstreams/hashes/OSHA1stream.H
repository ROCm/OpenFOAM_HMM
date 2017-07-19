/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::OSHA1stream

Description
    An output stream for calculating SHA1 digests.

SourceFiles
    OSHA1stream.C

\*---------------------------------------------------------------------------*/

#ifndef OSHA1stream_H
#define OSHA1stream_H

#include "OSstream.H"
#include "SHA1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class osha1stream;
class OSHA1stream;

/*---------------------------------------------------------------------------*\
                        Class sha1streambuf Declaration
\*---------------------------------------------------------------------------*/

//- A streambuf class for calculating SHA1 digests
class sha1streambuf
:
    public std::streambuf
{
    // Private data

    //- This does all the work and has its own buffering
    SHA1 sha1_;

    friend class osha1stream;

public:

    // Constructors

        //- Construct null
        sha1streambuf()
        {}

    // Member Functions

      // Write

        //- Process unbuffered
        virtual std::streamsize xsputn(const char* str, std::streamsize n)
        {
            sha1_.append(str, n);
            return n;
        }
};


/*---------------------------------------------------------------------------*\
                         Class osha1stream Declaration
\*---------------------------------------------------------------------------*/

//- A basic output stream for calculating SHA1 digests
class osha1stream
:
    virtual public std::ios,
    public std::ostream
{
    // Private data

    sha1streambuf sbuf_;

public:

    // Constructors

        //- Construct null
        osha1stream()
        :
            std::ostream(&sbuf_)
        {}

    // Member Functions

      // Access

        //- This hides both signatures of std::basic_ios::rdbuf()
        sha1streambuf* rdbuf()
        {
            return &sbuf_;
        }

        //- Full access to the sha1
        SHA1& sha1()
        {
            return sbuf_.sha1_;
        }

};


/*---------------------------------------------------------------------------*\
                    Class OSHA1streamAllocator Declaration
\*---------------------------------------------------------------------------*/

//- Allocator for an osha1stream
class OSHA1streamAllocator
{
protected:

    // Member data

        //- The allocated stream pointer
        osha1stream* allocatedPtr_;


    // Constructors

        //- Construct null
        OSHA1streamAllocator()
        :
            allocatedPtr_(new osha1stream())
        {}


    //- Destructor
    ~OSHA1streamAllocator()
    {
        deallocate();
    }


    // Member Functions

        //- Delete the stream pointer
        void deallocate()
        {
            if (allocatedPtr_)
            {
                delete allocatedPtr_;
                allocatedPtr_ = nullptr;
            }
        }


public:

        //- Full access to the sha1
        SHA1& sha1()
        {
            return allocatedPtr_->sha1();
        }

};


/*---------------------------------------------------------------------------*\
                         Class OSHA1stream Declaration
\*---------------------------------------------------------------------------*/

//- The output stream for calculating SHA1 digests
class OSHA1stream
:
    public OSHA1streamAllocator,
    public OSstream
{

    // Private Member Functions

        //- Disallow default bitwise copy construct
        OSHA1stream(const OSHA1stream&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const OSHA1stream&) = delete;

public:

    // Constructors

        //- Construct with an empty digest
        OSHA1stream
        (
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            OSHA1streamAllocator(),
            OSstream(*allocatedPtr_, "OSHA1stream", format, version)
        {}


    //- Destructor
    ~OSHA1stream()
    {}


    // Member functions

      // Access

        //- Return SHA1::Digest for the data processed until now
        SHA1Digest digest()
        {
            return sha1().digest();
        }


      // Edit

        //- Clear the SHA1 calculation
        void reset()
        {
            sha1().clear();
        }

        //- Clear the SHA1 calculation
        //  \deprecated use reset instead (deprecated Jul 2017)
        void rewind()
        {
            sha1().clear();
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
