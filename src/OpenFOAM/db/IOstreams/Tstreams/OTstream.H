/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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
    Foam::OTstream

Description
    A simple output token stream that can be used to build token lists.
    Always UNCOMPRESSED.

Note
    Appending single characters to token list is fragile.

SourceFiles
    OTstream.C

\*---------------------------------------------------------------------------*/

#ifndef OTstream_H
#define OTstream_H

#include "token.H"
#include "Ostream.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class OTstream Declaration
\*---------------------------------------------------------------------------*/

class OTstream
:
    public Ostream,
    public DynamicList<token>
{
public:

    // Constructors

        //- Default construct, set stream status
        explicit OTstream(IOstreamOption streamOpt = IOstreamOption())
        :
            Ostream(streamOpt.format(), streamOpt.version()),
            DynamicList<token>()
        {
            setOpened();
            setGood();
        }

        //- Copy construct
        OTstream(const OTstream& os)
        :
            Ostream(os.format(), os.version()),
            DynamicList<token>(os.tokens())
        {
            setOpened();
            setGood();
        }

        //- Move construct
        OTstream(OTstream&& os)
        :
            Ostream(os.format(), os.version()),
            DynamicList<token>(std::move(os.tokens()))
        {
            setOpened();
            setGood();
        }


    //- Destructor
    ~OTstream() = default;


    // Member Functions

        //- The tokens
        const DynamicList<token>& tokens() const
        {
            return *this;
        }

        //- The tokens
        DynamicList<token>& tokens()
        {
            return *this;
        }


    // Write

        //- Write token to stream or otherwise handle it.
        //  \return false if the token type was not handled by this method
        virtual bool write(const token& tok);

        //- Write single character. Whitespace is suppressed.
        virtual Ostream& write(const char c);

        //- Write the word-characters of a character string.
        //  Sends as a single char, or as word.
        virtual Ostream& write(const char* str);

        //- Write word
        virtual Ostream& write(const word& str);

        //- Write string
        virtual Ostream& write(const string& str);

        //- Write std::string surrounded by quotes.
        //  Optional write without quotes.
        virtual Ostream& writeQuoted
        (
            const std::string& str,
            const bool quoted=true
        );

        //- Write int32_t as a label
        virtual Ostream& write(const int32_t val);

        //- Write int64_t as a label
        virtual Ostream& write(const int64_t val);

        //- Write floatScalar
        virtual Ostream& write(const floatScalar val);

        //- Write doubleScalar
        virtual Ostream& write(const doubleScalar val);

        //- Write binary block with 8-byte alignment.
        virtual Ostream& write(const char* data, std::streamsize count);

        //- Low-level raw binary output.
        virtual Ostream& writeRaw(const char* data, std::streamsize count);

        //- Begin marker for low-level raw binary output.
        //  The count indicates the number of bytes for subsequent
        //  writeRaw calls.
        virtual bool beginRawWrite(std::streamsize count);

        //- End marker for low-level raw binary output.
        virtual bool endRawWrite()
        {
            return true;
        }

        //- Add indentation characters
        virtual void indent()
        {}


    // Stream State Functions

        //- Get flags of output stream
        virtual ios_base::fmtflags flags() const
        {
            return ios_base::fmtflags(0);
        }

        //- Set flags of stream - ignored
        std::ios_base::fmtflags flags(const ios_base::fmtflags)
        {
            return ios_base::fmtflags(0);
        }

        //- Flush stream
        virtual void flush()
        {}

        //- Add newline and flush stream
        virtual void endl()
        {}

        //- Get the current padding character
        //  \return previous padding character
        virtual char fill() const
        {
            return 0;
        }

        //- Set padding character for formatted field up to field width
        virtual char fill(const char)
        {
            return 0;
        }

        //- Get width of output field
        virtual int width() const
        {
            return 0;
        }

        //- Set width of output field
        //  \return previous width
        virtual int width(const int)
        {
             return 0;
        }

        //- Get precision of output field
        virtual int precision() const
        {
             return 0;
        }

        //- Set precision of output field
        //  \return old precision
        virtual int precision(const int)
        {
             return 0;
        }


    // Other

        //- Reset the output buffer and rewind the stream
        void reset()
        {
            this->rewind();
        }

        //- Rewind the output stream
        virtual void rewind()
        {
            DynamicList<token>::clear();
            setOpened();
            setGood();
        }

        //- Print stream description to Ostream
        void print(Ostream& os) const;


    // Additional constructors and methods (as per v2012 and earlier)
    #ifdef Foam_IOstream_extras

        //- Construct empty with format, version
        explicit OTstream
        (
            IOstreamOption::streamFormat fmt,
            IOstreamOption::versionNumber ver = IOstreamOption::currentVersion
        )
        :
            OTstream(IOstreamOption(fmt, ver))
        {}

    #endif /* Foam_IOstream_extras */
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
