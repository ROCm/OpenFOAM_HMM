/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::ISstream

Description
    Generic input stream.

SourceFiles
    ISstreamI.H
    ISstream.C

\*---------------------------------------------------------------------------*/

#ifndef ISstream_H
#define ISstream_H

#include "Istream.H"
#include "fileName.H"
#include <iostream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class ISstream Declaration
\*---------------------------------------------------------------------------*/

class ISstream
:
    public Istream
{
    // Private data

        fileName name_;
        std::istream& is_;


    // Private Member Functions

        char nextValid();

        void readWordToken(token& t);

    // Private Member Functions


        //- Read a verbatim string (excluding block delimiters).
        //  The leading "#{" has been removed prior to calling,
        //  continues until the closing "#}" has been found.
        Istream& readVerbatim(string& str);

        //- Read a variable name starting with '$'.
        //  Handles both "$var" and "${var}" forms.
        Istream& readVariable(string& str);

        //- Disallow default bitwise assignment
        void operator=(const ISstream&) = delete;


public:

    // Constructors

        //- Construct as wrapper around std::istream
        inline ISstream
        (
            std::istream& is,
            const string& name,
            streamFormat format=ASCII,
            versionNumber version=currentVersion,
            compressionType compression=UNCOMPRESSED
        );


    //- Destructor
    virtual ~ISstream()
    {}


    // Member functions

        // Inquiry

            //- Return the name of the stream
            //  Useful for Fstream to return the filename
            virtual const fileName& name() const
            {
                return name_;
            }

            //- Return non-const access to the name of the stream
            //  Useful to alter the stream name
            virtual fileName& name()
            {
                return name_;
            }

            //- Return flags of output stream
            virtual ios_base::fmtflags flags() const;


        // Read functions

            //- Raw, low-level get character function.
            inline ISstream& get(char& c);

            //- Raw, low-level peek function.
            //  Does not remove the character from the stream.
            //  Returns the next character in the stream or EOF if the
            //  end of file is read.
            inline int peek();

            //- Raw, low-level getline into a string function.
            inline ISstream& getLine(string& str);

            //- Raw, low-level putback character function.
            inline ISstream& putback(const char c);

            //- Return next token from stream
            virtual Istream& read(token& t);

            //- Read a character
            virtual Istream& read(char& c);

            //- Read a word
            virtual Istream& read(word& str);

            //- Read a string (including enclosing double-quotes).
            //  Backslashes are retained, except when escaping double-quotes
            //  and an embedded newline character.
            virtual Istream& read(string& str);

            //- Read a label
            virtual Istream& read(label& val);

            //- Read a floatScalar
            virtual Istream& read(floatScalar& val);

            //- Read a doubleScalar
            virtual Istream& read(doubleScalar& val);

            //- Read binary block
            virtual Istream& read(char* buf, std::streamsize count);

            //- Rewind the stream so that it may be read again
            virtual void rewind();


        // Stream state functions

            //- Set stream flags
            virtual ios_base::fmtflags flags(const ios_base::fmtflags flags);


        // STL stream

            //- Access to underlying std::istream
            virtual std::istream& stdStream()
            {
                return is_;
            }

            //- Const access to underlying std::istream
            virtual const std::istream& stdStream() const
            {
                return is_;
            }


        // Print

            //- Print description of IOstream to Ostream
            virtual void print(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ISstreamI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
