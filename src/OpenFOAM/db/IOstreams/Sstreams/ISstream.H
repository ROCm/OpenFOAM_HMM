/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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
    Generic input stream using a standard (STL) stream.

SourceFiles
    ISstreamI.H
    ISstream.C

\*---------------------------------------------------------------------------*/

#ifndef ISstream_H
#define ISstream_H

#include "Istream.H"
#include "fileName.H"
#include <limits>
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
    // Private Data

        fileName name_;

        std::istream& is_;


    // Private Member Functions

        //- Get the next valid (non-whitespace) character,
        //- after skipping any C/C++ comments.
        char nextValid();

        //- No copy assignment
        void operator=(const ISstream&) = delete;


public:

    // Constructors

        //- Construct wrapper around std::istream, set stream status
        //  Default stream options (ASCII, uncompressed)
        inline ISstream
        (
            std::istream& is,
            const string& streamName,
            IOstreamOption streamOpt = IOstreamOption()
        );

        //- Construct wrapper around std::istream, set stream status
        ISstream
        (
            std::istream& is,
            const string& streamName,
            IOstreamOption::streamFormat fmt,
            IOstreamOption::versionNumber ver = IOstreamOption::currentVersion,
            IOstreamOption::compressionType cmp = IOstreamOption::UNCOMPRESSED
        )
        :
            ISstream(is, streamName, IOstreamOption(fmt, ver, cmp))
        {}


    //- Destructor
    virtual ~ISstream() = default;


    // Member Functions

    // Characteristics

        //- Return the name of the stream.
        //  Useful for Fstream to return the filename
        virtual const fileName& name() const
        {
            return name_;
        }

        //- Return stream name for modification
        virtual fileName& name()
        {
            return name_;
        }

        //- Return flags of output stream
        virtual ios_base::fmtflags flags() const;


    // Special-purpose Functions

        //- Discard until end of C-style comment '*/'
        //  \return False if stream exhausted before finding the comment end
        bool seekCommentEnd_Cstyle();

        //- Raw, low-level get into a string.
        //- Continues reading \b after an initial left-brace until it finds
        //- the matching closing right-brace.
        //  Tracks balanced pairs, trims out leading/trailing whitespace.
        //
        //  \return False if stream exhausted before finding closing brace
        bool continueReadUntilRightBrace
        (
            std::string& str,
            const bool stripComments = true
        );


        // Read Functions

            //- Raw, low-level get character function.
            inline ISstream& get(char& c);

            //- Raw, low-level peek function.
            //  Does not remove the character from the stream.
            //  Returns the next character in the stream or EOF if the
            //  end of file is read.
            inline int peek();

            //- Raw, low-level getline (until delimiter) into a string.
            inline ISstream& getLine(std::string& str, char delim = '\n');

            //- Low-level discard until delimiter
            //  \return the number of characters extracted
            inline std::streamsize getLine(std::nullptr_t, char delim = '\n');

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

            //- Low-level raw binary read
            virtual Istream& readRaw(char* data, std::streamsize count);

            //- Start of low-level raw binary read
            virtual bool beginRawRead();

            //- End of low-level raw binary read
            virtual bool endRawRead();

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

            //- Print stream description to Ostream
            virtual void print(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ISstreamI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
