/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

#ifndef Foam_ITstream_H
#define Foam_ITstream_H

#include "Istream.H"
#include "tokenList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class ITstream Declaration
\*---------------------------------------------------------------------------*/

class ITstream
:
    public Istream,
    public tokenList
{
    // Private Data

        //- Name associated with the stream
        fileName name_;

        //- Index of token currently being read
        label tokenIndex_;


    // Private Member Functions

        //- An ad hoc combination of reserve and setCapacity somewhat
        //- similar to DynamicList.
        //
        //  In lazy mode, increase list size if needed, but leave any
        //  excess capacity - works like reserve.
        //
        //  In non-lazy mode, set exact capacity
        void reserveCapacity(const label nElem, const bool lazy);


        //- Failsafe read-access to token at specified location
        //- or undefinedToken
        inline const token& peekAt(const label i) const
        {
            return
            (
                i >= 0 && i < tokenList::size()
                ? tokenList::operator[](i)
                : token::undefinedToken
            );
        }


public:

    // Constructors

        //- Copy construct
        ITstream(const ITstream& is);

        //- Move construct
        ITstream(ITstream&& is);

        //- Default construct. Empty stream, optionally with given name
        explicit ITstream
        (
            IOstreamOption streamOpt = IOstreamOption(),
            const string& name = "input"
        );

        //- Construct empty, optionally with given name
        explicit ITstream
        (
            const Foam::zero,
            const string& name = "input",
            IOstreamOption streamOpt = IOstreamOption()
        );

        //- Copy construct from tokens, optionally with given name
        explicit ITstream
        (
            const UList<token>& tokens,
            IOstreamOption streamOpt = IOstreamOption(),
            const string& name = "input"
        );

        //- Move construct from tokens, optionally with given name
        explicit ITstream
        (
            List<token>&& tokens,
            IOstreamOption streamOpt = IOstreamOption(),
            const string& name = "input"
        );

        //- Construct token list by parsing the input character sequence
        //  Uses static parse function internally.
        explicit ITstream
        (
            const UList<char>& input,
            IOstreamOption streamOpt = IOstreamOption(),
            const string& name = "input"
        );

        //- Construct token list by parsing the input string
        //  Uses static parse function internally.
        explicit ITstream
        (
            const std::string& input,
            IOstreamOption streamOpt = IOstreamOption(),
            const string& name = "input"
        );

        //- Construct token list by parsing the input character sequence
        //  Uses static parse function internally.
        explicit ITstream
        (
            const char* input,
            IOstreamOption streamOpt = IOstreamOption(),
            const string& name = "input"
        );


    // Additional constructors

        //- Copy construct from tokens, with given name
        ITstream
        (
            const string& name,
            const UList<token>& tokens,
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            ITstream(tokens, streamOpt, name)
        {}

        //- Move construct from tokens, with given name
        ITstream
        (
            const string& name,
            List<token>&& tokens,
            IOstreamOption streamOpt = IOstreamOption()
        )
        :
            ITstream(std::move(tokens), streamOpt, name)
        {}


    //- Destructor
    virtual ~ITstream() = default;


    // Static Functions

        //- Create token list by parsing the input character sequence
        //- until no good tokens remain.
        static tokenList parse
        (
            const UList<char>& input,
            IOstreamOption streamOpt = IOstreamOption()
        );

        //- Create token list by parsing the input string
        //- until no good tokens remain.
        static tokenList parse
        (
            const std::string& input,
            IOstreamOption streamOpt = IOstreamOption()
        );

        //- Create token list by parsing the input character sequence
        //- until no good tokens remain.
        static tokenList parse
        (
            const char* input,
            IOstreamOption streamOpt = IOstreamOption()
        );


    // Member Functions

    // Characteristics

        //- The name of the input token stream.
        virtual const fileName& name() const { return name_; }

        //- The name of the input token stream, for modification.
        virtual fileName& name() { return name_; }


    // Token Access

        //- True if putback token is in use
        bool hasPutback() const noexcept
        {
            return Istream::hasPutback();
        }

        //- Failsafe peek at the \b first token in the list.
        //  \return \c undefinedToken if the list is empty.
        const token& front() const { return peekAt(0); }

        //- Failsafe peek at the \b last token in the list.
        //  \return \c undefinedToken if the list is empty.
        const token& back() const { return peekAt(tokenList::size()-1); }

        //- Failsafe peek at what the next read would return,
        //- including handling of any putback
        //  \return \c undefinedToken if list is exhausted
        const token& peek() const;

        //- The current token index when reading, or the insertion point.
        label tokenIndex() const noexcept { return tokenIndex_; }

        //- Non-const access to the current token index
        label& tokenIndex() noexcept { return tokenIndex_; }

        //- Number of tokens remaining
        label nRemainingTokens() const noexcept
        {
            return size() - tokenIndex_;
        }

        //- Move tokenIndex to the specified position
        //  Using seek(0) is identical to rewind.
        //  Using seek(-1) moves to the end.
        void seek(label pos);

        //- Move tokenIndex relative to the current position.
        //  Will not overrun the beginning or end positions.
        //
        //  Use skip(2) to move forward two tokens.
        //  Use skip(-2) to move backward two tokens.
        void skip(label n = 1);


    // Token list modification

        //- Copy append a token at the current tokenIndex,
        //- incrementing the index.
        void push_back(const token& t, const bool lazy);

        //- Move append a token at the current tokenIndex,
        //- incrementing the index.
        void push_back(token&& t, const bool lazy);

        //- Copy append a list of tokens at the current tokenIndex,
        //- incrementing the index.
        //
        //  \param newTokens the list of tokens to copy append
        //  \param lazy leaves any excess capacity for further appends.
        //      The caller will be responsible for resizing later.
        void push_back(const UList<token>& newTokens, const bool lazy);

        //- Move append a list of tokens at the current tokenIndex,
        //- incrementing the index.
        //
        //  \param newTokens the list of tokens to move append
        //  \param lazy leaves any excess capacity for further appends.
        //      The caller will be responsible for resizing later.
        void push_back(List<token>&& newTokens, const bool lazy);


    // Stream State Functions

        //- Get stream flags - always 0
        virtual ios_base::fmtflags flags() const
        {
            return ios_base::fmtflags(0);
        }

        //- Set flags of stream - ignored
        ios_base::fmtflags flags(const ios_base::fmtflags)
        {
            return ios_base::fmtflags(0);
        }


        // Read Functions

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

            //- Read a float
            virtual Istream& read(float&);

            //- Read a double
            virtual Istream& read(double&);

            //- Read binary block
            //  \note Not implemented
            virtual Istream& read(char* data, std::streamsize);

            //- Low-level raw binary read
            //  \note Not implemented
            virtual Istream& readRaw(char* data, std::streamsize count);

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
            virtual void rewind();


    // Output

        //- Print stream description to Ostream
        void print(Ostream& os) const;

        //- Concatenate tokens into a space-separated std::string.
        //- The resulting string may contain quote characters.
        std::string toString() const;


    // Member Operators

        //- Use operator() from Istream
        using Istream::operator();

        //- Copy assignment, with rewind()
        void operator=(const ITstream& is);

        //- Copy assignment of tokens, with rewind()
        void operator=(const UList<token>& toks);

        //- Move assignment of tokens, with rewind()
        void operator=(List<token>&& toks);


    // Additional constructors and methods (as per v2012 and earlier)
    #ifdef Foam_IOstream_extras

        //- Construct from components, copying the tokens
        ITstream
        (
            const string& name,
            const UList<token>& tokens,
            IOstreamOption::streamFormat fmt
        )
        :
            ITstream(name, tokens, IOstreamOption(fmt))
        {}

        //- Construct from components, transferring the tokens
        ITstream
        (
            const string& name,
            List<token>&& tokens,
            IOstreamOption::streamFormat fmt
        )
        :
            ITstream(name, std::move(tokens), IOstreamOption(fmt))
        {}

    #endif /* Foam_IOstream_extras */


    // Housekeeping

        //- Same as front()
        const token& peekFirst() const { return front(); }

        //- Copy append a token at the current tokenIndex,
        //- incrementing the index.
        void append(const token& t, const bool lazy)
        {
            this->push_back(std::move(t), lazy);
        }

        //- Move append a token at the current tokenIndex,
        //- incrementing the index.
        void append(token&& t, const bool lazy)
        {
            this->push_back(std::move(t), lazy);
        }

        //- Copy append a list of tokens at the current tokenIndex,
        //- incrementing the index.
        //
        //  \param newTokens the list of tokens to copy append
        //  \param lazy leaves any excess capacity for further appends.
        //      The caller will be responsible for resizing later.
        void append(const UList<token>& newTokens, const bool lazy)
        {
            this->push_back(newTokens, lazy);
        }

        //- Move append a list of tokens at the current tokenIndex,
        //- incrementing the index.
        //
        //  \param newTokens the list of tokens to move append
        //  \param lazy leaves any excess capacity for further appends.
        //      The caller will be responsible for resizing later.
        void append(List<token>&& newTokens, const bool lazy)
        {
            this->push_back(std::move(newTokens), lazy);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
