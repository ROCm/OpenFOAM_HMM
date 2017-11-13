/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
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

#include "UOPstream.H"
#include "int.H"
#include "token.H"
#include <cctype>

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline void Foam::UOPstream::prepareBuffer
(
    const size_t count,
    const size_t align
)
{
    if (!count)
    {
        return;
    }

    // The current output position
    label pos = sendBuf_.size();

    if (align > 1)
    {
        // Align output position. Pads sendBuf_.size() - oldPos characters.
        pos = align + ((pos - 1) & ~(align - 1));
    }

    // Extend buffer (as required)
    sendBuf_.reserve(max(1000, label(pos + count)));

    // Move to the aligned output position
    sendBuf_.setSize(pos);
}


template<class T>
inline void Foam::UOPstream::writeToBuffer(const T& val)
{
    writeToBuffer(&val, sizeof(T), sizeof(T));
}


inline void Foam::UOPstream::writeToBuffer(const char& c)
{
    if (!sendBuf_.capacity())
    {
        sendBuf_.setCapacity(1000);
    }
    sendBuf_.append(c);
}


inline void Foam::UOPstream::writeToBuffer
(
    const void* data,
    const size_t count,
    const size_t align
)
{
    if (!count)
    {
        return;
    }

    prepareBuffer(count, align);

    // The aligned output position
    const label pos = sendBuf_.size();

    // Extend the addressable range for direct pointer access
    sendBuf_.setSize(pos + count);

    char* const __restrict__ buf = (sendBuf_.begin() + pos);
    const char* const __restrict__ input = reinterpret_cast<const char*>(data);

    for (size_t i = 0; i < count; ++i)
    {
        buf[i] = input[i];
    }
}


inline void Foam::UOPstream::writeStringToBuffer(const std::string& str)
{
    const size_t len = str.size();
    writeToBuffer(len);
    writeToBuffer(str.data(), len, 1);
}


// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::UOPstream::UOPstream
(
    const commsTypes commsType,
    const int toProcNo,
    DynamicList<char>& sendBuf,
    const int tag,
    const label comm,
    const bool sendAtDestruct,
    streamFormat format,
    versionNumber version
)
:
    UPstream(commsType),
    Ostream(format, version),
    toProcNo_(toProcNo),
    sendBuf_(sendBuf),
    tag_(tag),
    comm_(comm),
    sendAtDestruct_(sendAtDestruct)
{
    setOpened();
    setGood();
}


Foam::UOPstream::UOPstream(const int toProcNo, PstreamBuffers& buffers)
:
    UPstream(buffers.commsType_),
    Ostream(buffers.format_, buffers.version_),
    toProcNo_(toProcNo),
    sendBuf_(buffers.sendBuf_[toProcNo]),
    tag_(buffers.tag_),
    comm_(buffers.comm_),
    sendAtDestruct_(buffers.commsType_ != UPstream::commsTypes::nonBlocking)
{
    setOpened();
    setGood();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::UOPstream::~UOPstream()
{
    if (sendAtDestruct_)
    {
        if
        (
            !UOPstream::write
            (
                commsType_,
                toProcNo_,
                sendBuf_.begin(),
                sendBuf_.size(),
                tag_,
                comm_
            )
        )
        {
            FatalErrorInFunction
                << "Failed sending outgoing message of size " << sendBuf_.size()
                << " to processor " << toProcNo_
                << Foam::abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::UOPstream::write(const token& tok)
{
    // Direct token handling only for some types

    switch (tok.type())
    {
        case token::tokenType::FLAG :
        {
            writeToBuffer(char(token::tokenType::FLAG));
            writeToBuffer(char(tok.flagToken()));

            return true;
        }

        case token::tokenType::VERBATIMSTRING :
        {
            writeToBuffer(char(token::tokenType::VERBATIMSTRING));
            write(tok.stringToken());

            return true;
        }

        case token::tokenType::VARIABLE :
        {
            writeToBuffer(char(token::tokenType::VARIABLE));
            write(tok.stringToken());

            return true;
        }

        default:
            break;
    }

    return false;
}


Foam::Ostream& Foam::UOPstream::write(const char c)
{
    if (!isspace(c))
    {
        writeToBuffer(c);
    }

    return *this;
}


Foam::Ostream& Foam::UOPstream::write(const char* str)
{
    const word nonWhiteChars(string::validate<word>(str));

    if (nonWhiteChars.size() == 1)
    {
        return write(nonWhiteChars[0]);
    }
    else if (nonWhiteChars.size())
    {
        return write(nonWhiteChars);
    }
    else
    {
        return *this;
    }
}


Foam::Ostream& Foam::UOPstream::write(const word& str)
{
    writeToBuffer(char(token::tokenType::WORD));
    writeStringToBuffer(str);

    return *this;
}


Foam::Ostream& Foam::UOPstream::write(const string& str)
{
    writeToBuffer(char(token::tokenType::STRING));
    writeStringToBuffer(str);

    return *this;
}


Foam::Ostream& Foam::UOPstream::writeQuoted
(
    const std::string& str,
    const bool quoted
)
{
    if (quoted)
    {
        writeToBuffer(char(token::tokenType::STRING));
    }
    else
    {
        writeToBuffer(char(token::tokenType::WORD));
    }
    writeStringToBuffer(str);

    return *this;
}


Foam::Ostream& Foam::UOPstream::write(const int32_t val)
{
    writeToBuffer(char(token::tokenType::LABEL));
    writeToBuffer(val);
    return *this;
}


Foam::Ostream& Foam::UOPstream::write(const int64_t val)
{
    writeToBuffer(char(token::tokenType::LABEL));
    writeToBuffer(val);
    return *this;
}


Foam::Ostream& Foam::UOPstream::write(const floatScalar val)
{
    writeToBuffer(char(token::tokenType::FLOAT_SCALAR));
    writeToBuffer(val);
    return *this;
}


Foam::Ostream& Foam::UOPstream::write(const doubleScalar val)
{
    writeToBuffer(char(token::tokenType::DOUBLE_SCALAR));
    writeToBuffer(val);
    return *this;
}


Foam::Ostream& Foam::UOPstream::write
(
    const char* data,
    const std::streamsize count
)
{
    if (format() != BINARY)
    {
        FatalErrorInFunction
            << "stream format not binary"
            << Foam::abort(FatalError);
    }

    writeToBuffer(data, count, 8);

    return *this;
}


Foam::Ostream& Foam::UOPstream::beginRaw
(
    const std::streamsize count
)
{
    if (format() != BINARY)
    {
        FatalErrorInFunction
            << "stream format not binary"
            << Foam::abort(FatalError);
    }

    // Alignment = 8, as per write(const char*, streamsize)
    prepareBuffer(count, 8);

    return *this;
}


Foam::Ostream& Foam::UOPstream::writeRaw
(
    const char* data,
    const std::streamsize count
)
{
    // No check for format() == BINARY since this is either done in the
    // beginRaw() method, or the caller knows what they are doing.

    // Previously aligned and sizes reserved via beginRaw()
    writeToBuffer(data, count, 1);

    return *this;
}


void Foam::UOPstream::print(Ostream& os) const
{
    os  << "Writing from processor " << toProcNo_
        << " to processor " << myProcNo() << " in communicator " << comm_
        << " and tag " << tag_ << Foam::endl;
}


// ************************************************************************* //
