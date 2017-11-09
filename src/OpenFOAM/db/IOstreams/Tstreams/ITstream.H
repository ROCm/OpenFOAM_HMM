/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::ITstream

Description
    An input stream of tokens.

SourceFiles
    ITstream.C

\*---------------------------------------------------------------------------*/

#ifndef ITstream_H
#define ITstream_H

#include "Istream.H"
#include "tokenList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration
class ISstream;

/*---------------------------------------------------------------------------*\
                           Class ITstream Declaration
\*---------------------------------------------------------------------------*/

class ITstream
:
    public Istream,
    public tokenList
{
    // Private data

        //- Name of ITstream
        fileName name_;

        //- Index of token currently being read
        label tokenIndex_;


    // Private Member Functions

        //- Convert input sequence into a list of tokens.
        //  \return the number of tokens in the resulting list.
        static label parseStream(ISstream& input, tokenList& tokens);


public:

    // Constructors

        //- Construct from components
        ITstream
        (
            const string& name,
            const UList<token>& tokens,
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            Istream(format, version),
            tokenList(tokens),
            name_(name),
            tokenIndex_(0)
        {
            setOpened();
            setGood();
        }


        //- Construct from components, transferring the tokens
        ITstream
        (
            const string& name,
            const Xfer<List<token>>& tokens,
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        )
        :
            Istream(format, version),
            tokenList(tokens),
            name_(name),
            tokenIndex_(0)
        {
            setOpened();
            setGood();
        }


        //- Construct token list by parsing the input character sequence
        //  Uses UIListStream internally.
        ITstream
        (
            const string& name,
            const UList<char>& input,
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        );


        //- Construct token list by parsing the input string
        //  Uses UIListStream internally.
        ITstream
        (
            const string& name,
            const std::string& input,
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        );


        //- Construct token list by parsing the input character sequence
        //  Uses UIListStream internally.
        ITstream
        (
            const string& name,
            const char* input,
            streamFormat format=ASCII,
            versionNumber version=currentVersion
        );


        //- Construct as copy
        ITstream(const ITstream& its)
        :
            Istream(ASCII, currentVersion),
            tokenList(its),
            name_(its.name_),
            tokenIndex_(0)
        {
            setOpened();
            setGood();
        }


        //- Destructor
        virtual ~ITstream()
        {}


    // Static Functions

        //- Create token list by parsing the input character sequence until
        //- no good tokens remain.
        static tokenList parse
        (
            const UList<char>& input,
            streamFormat format=ASCII
        );

        //- Create token list by parsing the input string until
        //- no good tokens remain.
        static tokenList parse
        (
            const std::string& input,
            streamFormat format=ASCII
        );

        //- Create token list by parsing the input character sequence until
        //- no good tokens remain.
        static tokenList parse
        (
            const char* input,
            streamFormat format=ASCII
        );


    // Member functions

        // Inquiry

            //- Return the name of the stream
            virtual const fileName& name() const
            {
                return name_;
            }

            //- Return non-const access to the name of the stream
            virtual fileName& name()
            {
                return name_;
            }

            //- Return the current token index
            label tokenIndex() const
            {
                return tokenIndex_;
            }

            //- Return non-const access to the current token index
            label& tokenIndex()
            {
                return tokenIndex_;
            }

            //- Return the number of remaining tokens
            label nRemainingTokens() const
            {
                return size() - tokenIndex_;
            }

            //- Return flags of output stream
            ios_base::fmtflags flags() const
            {
                return ios_base::fmtflags(0);
            }


        // Read functions

            //- Return next token from stream
            virtual Istream& read(token& tok);

            //- Read a character
            virtual Istream& read(char&);

            //- Read a word
            virtual Istream& read(word&);

            // Read a string (including enclosing double-quotes)
            virtual Istream& read(string&);

            //- Read a label
            virtual Istream& read(label&);

            //- Read a floatScalar
            virtual Istream& read(floatScalar&);

            //- Read a doubleScalar
            virtual Istream& read(doubleScalar&);

            //- Read binary block
            virtual Istream& read(char*, std::streamsize);

            //- Rewind the stream so that it may be read again
            virtual void rewind();


        // Edit

            //- Set flags of stream
            ios_base::fmtflags flags(const ios_base::fmtflags)
            {
                return ios_base::fmtflags(0);
            }


        // Print

            //- Print description of stream to Ostream
            void print(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
