/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Foam::dummyISstream

Description
    Dummy input stream, which can be used as a placeholder for interfaces
    taking an Istream or ISstream. Aborts at any attempt to read from it.

Note
    The inheritance from an empty IStringStream is solely for convenience
    of implementation and should not be relied upon.

SourceFiles
    dummyISstream.H

\*---------------------------------------------------------------------------*/

#ifndef dummyISstream_H
#define dummyISstream_H

#include "StringStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class dummyISstream Declaration
\*---------------------------------------------------------------------------*/

class dummyISstream
:
    public IStringStream
{
public:

    // Generated Methods

        //- Default construct
        dummyISstream() = default;

        //- Destructor
        virtual ~dummyISstream() = default;


    // Member Functions

        // Read Functions

            //- Return next token from stream
            virtual Istream& read(token&)
            {
                NotImplemented;
                return *this;
            }

            //- Read a character
            virtual Istream& read(char&)
            {
                NotImplemented;
                return *this;
            }

            //- Read a word
            virtual Istream& read(word&)
            {
                NotImplemented;
                return *this;
            }

            // Read a string (including enclosing double-quotes)
            virtual Istream& read(string&)
            {
                NotImplemented;
                return *this;
            }

            //- Read a label
            virtual Istream& read(label&)
            {
                NotImplemented;
                return *this;
            }

            //- Read a floatScalar
            virtual Istream& read(floatScalar&)
            {
                NotImplemented;
                return *this;
            }

            //- Read a doubleScalar
            virtual Istream& read(doubleScalar&)
            {
                NotImplemented;
                return *this;
            }

            //- Read binary block
            virtual Istream& read(char*, std::streamsize)
            {
                NotImplemented;
                return *this;
            }

            //- Low-level raw binary read
            virtual Istream& readRaw(char*, std::streamsize)
            {
                NotImplemented;
                return *this;
            }

            //- Start of low-level raw binary read
            virtual bool beginRawRead()
            {
                return false;
            }

            //- End of low-level raw binary read
            virtual bool endRawRead()
            {
                return false;
            }

            //- Rewind the stream so that it may be read again
            virtual void rewind()
            {}

            //- Return flags of stream
            virtual ios_base::fmtflags flags() const
            {
                return ios_base::fmtflags(0);
            }

            //- Set flags of stream
            virtual ios_base::fmtflags flags(const ios_base::fmtflags)
            {
                return ios_base::fmtflags(0);
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
