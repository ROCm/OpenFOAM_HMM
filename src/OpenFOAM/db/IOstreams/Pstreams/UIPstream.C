/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "UIPstream.H"
#include "int.H"
#include "token.H"
#include <cctype>

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Convert a single character to a word with length 1
inline static Foam::word charToWord(char c)
{
    return Foam::word(std::string(1, c), false);
}


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


// Return the position with word boundary alignment
inline static label byteAlign(const label pos, const size_t align)
{
    return
    (
        (align > 1)
      ? (align + ((pos - 1) & ~(align - 1)))
      : pos
    );
}

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline void Foam::UIPstream::checkEof()
{
    if (recvBufPos_ == messageSize_)
    {
        setEof();
    }
}


inline void Foam::UIPstream::prepareBuffer(const size_t align)
{
    recvBufPos_ = byteAlign(recvBufPos_, align);
}


template<class T>
inline void Foam::UIPstream::readFromBuffer(T& val)
{
    prepareBuffer(sizeof(T));

    val = reinterpret_cast<T&>(recvBuf_[recvBufPos_]);
    recvBufPos_ += sizeof(T);
    checkEof();
}


inline void Foam::UIPstream::readFromBuffer
(
    void* data,
    const size_t count
)
{
    const char* const __restrict__ buf = &recvBuf_[recvBufPos_];
    char* const __restrict__ output = reinterpret_cast<char*>(data);

    for (size_t i = 0; i < count; ++i)
    {
        output[i] = buf[i];
    }

    recvBufPos_ += count;
    checkEof();
}


inline Foam::Istream& Foam::UIPstream::readString(std::string& str)
{
    // Use std::string::assign() to copy content, including '\0'.
    // Stripping (when desired) is the responsibility of the sending side.

    size_t len;
    readFromBuffer(len);

    if (len)
    {
        str.assign(&recvBuf_[recvBufPos_], len);
        recvBufPos_ += len;
        checkEof();
    }
    else
    {
        str.clear();
    }

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
                << " clearing receive buffer of size "
                << recvBuf_.size()
                << " messageSize_:" << messageSize_ << endl;
        }
        recvBuf_.clearStorage();
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
    t.lineNumber(this->lineNumber());

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
        case token::PLUS :
        case token::MINUS :
        case token::MULTIPLY :
        case token::DIVIDE :
        {
            t = token::punctuationToken(c);
            return *this;
        }

        // The word-variants
        case token::tokenType::WORD :
        case token::tokenType::DIRECTIVE :
        {
            word val;
            if (readString(val))
            {
                if (token::compound::isCompound(val))
                {
                    t = token::compound::New(val, *this).ptr();
                }
                else
                {
                    t = std::move(val);
                    t.setType(token::tokenType(c));
                }
            }
            else
            {
                t.setBad();
            }
            return *this;
        }

        // The string-variants
        case token::tokenType::STRING :
        case token::tokenType::EXPRESSION :
        case token::tokenType::VARIABLE :
        case token::tokenType::VERBATIM :
        {
            string val;
            if (readString(val))
            {
                t = std::move(val);
                t.setType(token::tokenType(c));
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

        // Float
        case token::tokenType::FLOAT :
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

        // Double
        case token::tokenType::DOUBLE :
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
                t = charToWord(c);
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
    c = recvBuf_[recvBufPos_];
    ++recvBufPos_;
    checkEof();
    return *this;
}


Foam::Istream& Foam::UIPstream::read(word& str)
{
    return readString(str);
}


Foam::Istream& Foam::UIPstream::read(string& str)
{
    return readString(str);
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
    if (count)
    {
        // For count == 0, a no-op
        // - see UOPstream::write(const char*, streamsize)
        beginRawRead();
        readRaw(data, count);
        endRawRead();
    }

    return *this;
}


Foam::Istream& Foam::UIPstream::readRaw(char* data, std::streamsize count)
{
    // No check for format() == BINARY since this is either done in the
    // beginRawRead() method, or the caller knows what they are doing.

    // Any alignment must have been done prior to this call
    readFromBuffer(data, count);
    return *this;
}


bool Foam::UIPstream::beginRawRead()
{
    if (format() != BINARY)
    {
        FatalErrorInFunction
            << "stream format not binary"
            << Foam::abort(FatalError);
    }

    // Align on word boundary (64-bit)
    // - as per read(const char*, streamsize)
    // The check for zero-size will have been done by the caller
    prepareBuffer(8);

    return true;
}


void Foam::UIPstream::rewind()
{
    recvBufPos_ = 0;
}


void Foam::UIPstream::print(Ostream& os) const
{
    os  << "Reading from processor " << fromProcNo_
        << " using communicator " << comm_
        <<  " and tag " << tag_
        << Foam::endl;
}


// ************************************************************************* //
