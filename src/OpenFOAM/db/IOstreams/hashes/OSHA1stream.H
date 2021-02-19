/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#ifndef OSHA1stream_H
#define OSHA1stream_H

#include "OSstream.H"
#include "SHA1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class osha1stream Declaration
\*---------------------------------------------------------------------------*/

//- A basic output stream for calculating SHA1 digests
class osha1stream
:
    virtual public std::ios,
    public std::ostream
{
    //- A streambuf class for calculating SHA1 digests
    class sha1buf
    :
        public std::streambuf
    {
        //- This does all the work and has its own buffering
        SHA1 sha1_;

    protected:

        //- Put sequence of characters
        virtual std::streamsize xsputn(const char* s, std::streamsize n)
        {
            sha1_.append(s, n);
            return n;
        }

    public:

        //- Default construct
        sha1buf() = default;

        //- Full access to the sha1
        inline SHA1& sha1()
        {
            return sha1_;
        }
    };


    // Private Data

        //- Reference to the underlying buffer
        sha1buf buf_;

public:

    // Constructors

        //- Default construct
        osha1stream()
        :
            std::ostream(&buf_)
        {}


    // Member Functions

        //- This hides both signatures of std::basic_ios::rdbuf()
        sha1buf* rdbuf()
        {
            return &buf_;
        }

        //- Full access to the sha1
        SHA1& sha1()
        {
            return buf_.sha1();
        }
};


namespace Detail
{

/*---------------------------------------------------------------------------*\
                Class Detail::OSHA1streamAllocator Declaration
\*---------------------------------------------------------------------------*/

//- Allocator for an osha1stream
class OSHA1streamAllocator
{
protected:

    // Protected Data

        typedef osha1stream stream_type;

        //- The output stream
        stream_type stream_;


    // Constructors

        //- Default construct
        OSHA1streamAllocator() = default;


public:

    // Member Functions

        //- Full access to the sha1
        SHA1& sha1()
        {
            return stream_.sha1();
        }


        //- Return SHA1::Digest for the data processed until now
        SHA1Digest digest()
        {
            return stream_.sha1().digest();
        }


        //- Clear the SHA1 calculation
        void reset()
        {
            return stream_.sha1().clear();
        }
};

} // End namespace Detail


/*---------------------------------------------------------------------------*\
                         Class OSHA1stream Declaration
\*---------------------------------------------------------------------------*/

//- The output stream for calculating SHA1 digests
class OSHA1stream
:
    public Detail::OSHA1streamAllocator,
    public OSstream
{
    typedef Detail::OSHA1streamAllocator allocator_type;

    // Private Member Functions

        //- No copy construct
        OSHA1stream(const OSHA1stream&) = delete;

        //- No copy assignment
        void operator=(const OSHA1stream&) = delete;

public:

    // Constructors

        //- Construct with an empty digest
        explicit OSHA1stream
        (
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            allocator_type(),
            OSstream(stream_, "sha1", streamOpt.format(), streamOpt.version())
        {}


    // Write Functions

        //- Add (unquoted) string contents.
        //  Ensures that SHA1 of C-string or C++-string content are identical.
        virtual Ostream& write(const string& str)
        {
            return writeQuoted(str, false); // Unquoted!
        }


    // Housekeeping

        //- Deprecated(2017-07) clear the SHA1 calculation
        //  \deprecated(2017-07) - use reset() method
        void rewind()
        {
            sha1().clear();
        }


    // Additional constructors and methods (as per v2012 and earlier)
    #ifdef Foam_IOstream_extras

        //- Construct empty
        explicit OSHA1stream
        (
            IOstreamOption::streamFormat fmt,
            IOstreamOption::versionNumber ver = IOstreamOption::currentVersion
        )
        :
            OSHA1stream(IOstreamOption(fmt, ver))
        {}

    #endif /* Foam_IOstream_extras */
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
