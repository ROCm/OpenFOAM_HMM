/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "UIPstream.H"
#include "int.H"
#include "token.H"
#include <cctype>

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
// Adjust stream format based on the flagMask
inline static void processFlags(Istream& is, int flagMask)
{
    if ((flagMask & token::ASCII))
    {
        is.format(IOstream::ASCII);
    }
    else if ((flagMask & token::BINARY))
    {
        is.format(IOstream::BINARY);
    }
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline void Foam::UIPstream::checkEof()
{
    if (externalBufPosition_ == messageSize_)
    {
        setEof();
    }
}


template<class T>
inline void Foam::UIPstream::readFromBuffer(T& val)
{
    const size_t align = sizeof(T);
    externalBufPosition_ = align + ((externalBufPosition_ - 1) & ~(align - 1));

    val = reinterpret_cast<T&>(externalBuf_[externalBufPosition_]);
    externalBufPosition_ += sizeof(T);
    checkEof();
}


inline void Foam::UIPstream::readFromBuffer
(
    void* data,
    const size_t count,
    const size_t align
)
{
    if (align > 1)
    {
        externalBufPosition_ =
            align
          + ((externalBufPosition_ - 1) & ~(align - 1));
    }

    const char* const __restrict__ buf = &externalBuf_[externalBufPosition_];
    char* const __restrict__ output = reinterpret_cast<char*>(data);

    for (size_t i = 0; i < count; ++i)
    {
        output[i] = buf[i];
    }

    externalBufPosition_ += count;
    checkEof();
}


inline Foam::Istream& Foam::UIPstream::readStringFromBuffer(std::string& str)
{
    // Use std::string::assign() to copy content, including '\0'.
    // Stripping (when desired) is the responsibility of the sending side.

    size_t len;
    readFromBuffer(len);

    if (len == 0)
    {
        str.clear();
    }
    else
    {
        str.assign(&externalBuf_[externalBufPosition_], len);
    }

    externalBufPosition_ += len;
    checkEof();

    return *this;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::UIPstream::~UIPstream()
{
    if (clearAtEnd_ && eof())
    {
        if (debug)
        {
            Pout<< "UIPstream::~UIPstream() : tag:" << tag_
                << " fromProcNo:" << fromProcNo_
                << " clearing externalBuf_ of size "
                << externalBuf_.size()
                << " messageSize_:" << messageSize_ << endl;
        }
        externalBuf_.clearStorage();
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Istream& Foam::UIPstream::read(token& t)
{
    // Return the put back token if it exists
    // - with additional handling for special stream flags
    if (Istream::getBack(t))
    {
        if (t.isFlag())
        {
            processFlags(*this, t.flagToken());
        }
        else
        {
            return *this;
        }
    }

    // Read character, return on error
    // - with additional handling for special stream flags

    char c;
    do
    {
        if (!read(c))
        {
            t.setBad();   // Error
            return *this;
        }

        if (c == token::FLAG)
        {
            char flagVal;

            if (read(flagVal))
            {
                processFlags(*this, flagVal);
            }
            else
            {
                t.setBad();   // Error
                return *this;
            }
        }
    }
    while (c == token::FLAG);


    // Set the line number of this token to the current stream line number
    t.lineNumber() = lineNumber();

    // Analyse input starting with this character.
    switch (c)
    {
        // Punctuation
        case token::END_STATEMENT :
        case token::BEGIN_LIST :
        case token::END_LIST :
        case token::BEGIN_SQR :
        case token::END_SQR :
        case token::BEGIN_BLOCK :
        case token::END_BLOCK :
        case token::COLON :
        case token::COMMA :
        case token::ASSIGN :
        case token::ADD :
        case token::SUBTRACT :
        case token::MULTIPLY :
        case token::DIVIDE :
        {
            t = token::punctuationToken(c);
            return *this;
        }

        // Word
        case token::tokenType::WORD :
        {
            word val;
            if (read(val))
            {
                if (token::compound::isCompound(val))
                {
                    t = token::compound::New(val, *this).ptr();
                }
                else
                {
                    t = std::move(val);
                }
            }
            else
            {
                t.setBad();
            }
            return *this;
        }

        // String
        case token::tokenType::VERBATIMSTRING :
        {
            // Recurse to read actual string
            read(t);
            t.setType(token::tokenType::VERBATIMSTRING);
            return *this;
        }
        case token::tokenType::VARIABLE :
        {
            // Recurse to read actual string
            read(t);
            t.setType(token::tokenType::VARIABLE);
            return *this;
        }
        case token::tokenType::STRING :
        {
            string val;
            if (read(val))
            {
                t = std::move(val);
            }
            else
            {
                t.setBad();
            }
            return *this;
        }

        // Label
        case token::tokenType::LABEL :
        {
            label val;
            if (read(val))
            {
                t = val;
            }
            else
            {
                t.setBad();
            }
            return *this;
        }

        // floatScalar
        case token::tokenType::FLOAT_SCALAR :
        {
            floatScalar val;
            if (read(val))
            {
                t = val;
            }
            else
            {
                t.setBad();
            }
            return *this;
        }

        // doubleScalar
        case token::tokenType::DOUBLE_SCALAR :
        {
            doubleScalar val;
            if (read(val))
            {
                t = val;
            }
            else
            {
                t.setBad();
            }
            return *this;
        }

        // Character (returned as a single character word) or error
        default:
        {
            if (isalpha(c))
            {
                t = word(c);
                return *this;
            }

            setBad();
            t.setBad();

            return *this;
        }
    }
}


Foam::Istream& Foam::UIPstream::read(char& c)
{
    c = externalBuf_[externalBufPosition_];
    ++externalBufPosition_;
    checkEof();
    return *this;
}


Foam::Istream& Foam::UIPstream::read(word& str)
{
    return readStringFromBuffer(str);
}


Foam::Istream& Foam::UIPstream::read(string& str)
{
    return readStringFromBuffer(str);
}


Foam::Istream& Foam::UIPstream::read(label& val)
{
    readFromBuffer(val);
    return *this;
}


Foam::Istream& Foam::UIPstream::read(floatScalar& val)
{
    readFromBuffer(val);
    return *this;
}


Foam::Istream& Foam::UIPstream::read(doubleScalar& val)
{
    readFromBuffer(val);
    return *this;
}


Foam::Istream& Foam::UIPstream::read(char* data, std::streamsize count)
{
    if (format() != BINARY)
    {
        FatalErrorInFunction
            << "stream format not binary"
            << Foam::abort(FatalError);
    }

    readFromBuffer(data, count, 8);
    return *this;
}


void Foam::UIPstream::rewind()
{
    externalBufPosition_ = 0;
}


void Foam::UIPstream::print(Ostream& os) const
{
    os  << "Reading from processor " << fromProcNo_
        << " using communicator " << comm_
        <<  " and tag " << tag_
        << Foam::endl;
}


// ************************************************************************* //
