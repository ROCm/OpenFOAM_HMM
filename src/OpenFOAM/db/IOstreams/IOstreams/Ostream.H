/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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
    Foam::Ostream

Description
    An Ostream is an abstract base class for all output systems
    (streams, files, token lists, etc).

SourceFiles
    Ostream.C

\*---------------------------------------------------------------------------*/

#ifndef Ostream_H
#define Ostream_H

#include "IOstream.H"
#include "keyType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class token;

/*---------------------------------------------------------------------------*\
                           Class Ostream Declaration
\*---------------------------------------------------------------------------*/

class Ostream
:
    public IOstream
{
protected:

    // Protected Data

        //- Indentation of the entry from the start of the keyword
        static constexpr const unsigned short entryIndentation_ = 16;

        //- Number of spaces per indent level
        unsigned short indentSize_ = 4;

        //- Current indent level
        unsigned short indentLevel_ = 0;


public:

    // Generated Methods

        //- Copy construct
        Ostream(const Ostream&) = default;

        //- Destructor
        virtual ~Ostream() = default;


    // Constructors

        //- Default construct (ASCII, uncompressed),
        //- construct with specified stream option
        explicit Ostream(IOstreamOption streamOpt = IOstreamOption())
        :
            IOstream(streamOpt)
        {}


        //- Construct with format, version (compression)
        explicit Ostream
        (
            IOstreamOption::streamFormat fmt,
            IOstreamOption::versionNumber ver = IOstreamOption::currentVersion,
            IOstreamOption::compressionType cmp = IOstreamOption::UNCOMPRESSED
        )
        :
            Ostream(IOstreamOption(fmt, ver, cmp))
        {}


    // Member Functions

        // Write Functions

            //- Write token to stream or otherwise handle it.
            //  \return false if the token type was not handled by this method
            virtual bool write(const token& tok) = 0;

            //- Write character
            virtual Ostream& write(const char c) = 0;

            //- Write character string
            virtual Ostream& write(const char* str) = 0;

            //- Write word
            virtual Ostream& write(const word& str) = 0;

            //- Write keyType
            //  A plain word is written unquoted.
            //  A regular expression is written as a quoted string.
            virtual Ostream& write(const keyType& kw);

            //- Write string
            virtual Ostream& write(const string& str) = 0;

            //- Write std::string surrounded by quotes.
            //  Optional write without quotes.
            virtual Ostream& writeQuoted
            (
                const std::string& str,
                const bool quoted=true
            ) = 0;

            //- Write int32_t
            virtual Ostream& write(const int32_t val) = 0;

            //- Write int64_t
            virtual Ostream& write(const int64_t val) = 0;

            //- Write floatScalar
            virtual Ostream& write(const floatScalar val) = 0;

            //- Write doubleScalar
            virtual Ostream& write(const doubleScalar val) = 0;

            //- Write binary block.
            virtual Ostream& write(const char* data, std::streamsize count) = 0;

            //- Low-level raw binary output.
            virtual Ostream& writeRaw
            (
                const char* data,
                std::streamsize count
            ) = 0;

            //- Emit begin marker for low-level raw binary output.
            //  The count indicates the number of bytes for subsequent
            //  writeRaw calls.
            virtual bool beginRawWrite(std::streamsize count) = 0;

            //- Emit end marker for low-level raw binary output.
            virtual bool endRawWrite() = 0;

            //- Add indentation characters
            virtual void indent() = 0;

            //- Return indent size (spaces per level)
            unsigned short indentSize() const noexcept
            {
                return indentSize_;
            }

            //- Change indent size (spaces per level), return old value
            unsigned short indentSize(unsigned short val) noexcept
            {
                auto old(indentSize_);
                indentSize_ = val;
                return old;
            }

            //- Return the indent level
            unsigned short indentLevel() const noexcept
            {
                return indentLevel_;
            }

            //- Change the indent level, return old value
            unsigned short indentLevel(unsigned short val) noexcept
            {
                auto old(indentLevel_);
                indentLevel_ = val;
                return old;
            }

            //- Increment the indent level
            void incrIndent() noexcept
            {
                ++indentLevel_;
            }

            //- Decrement the indent level
            void decrIndent();

            //- Write the keyword followed by an appropriate indentation
            virtual Ostream& writeKeyword(const keyType& kw);

            //- Write begin block group with the given name
            //  Increments indentation, adds newline.
            virtual Ostream& beginBlock(const keyType& kw);

            //- Write begin block group without a name
            //  Increments indentation, adds newline.
            virtual Ostream& beginBlock();

            //- Write end block group
            //  Decrements indentation, adds newline.
            virtual Ostream& endBlock();

            //- Write end entry (';') followed by newline.
            virtual Ostream& endEntry();

            //- Write a keyword/value entry.
            //  The following two are functionally equivalent:
            // \code
            // os.writeEntry(key, value);
            //
            // os.writeKeyword(key) << value << endEntry;
            // \endcode
            template<class T>
            Ostream& writeEntry(const keyType& key, const T& value)
            {
                writeKeyword(key) << value;
                return endEntry();
            }

            //- Write a keyword/value entry only when the two values differ.
            //  \param key the name of the entry
            //  \param value1 the reference value
            //  \param value2 the value to write if it differs from value1
            template<class T>
            Ostream& writeEntryIfDifferent
            (
                const word& key,
                const T& value1,
                const T& value2
            )
            {
                if (value1 != value2)
                {
                    writeEntry(key, value2);
                }

                return *this;
            }


        // Stream state functions

            //- Flush stream
            virtual void flush() = 0;

            //- Add newline and flush stream
            virtual void endl() = 0;

            //- Get padding character
            virtual char fill() const = 0;

            //- Set padding character for formatted field up to field width
            virtual char fill(const char fillch) = 0;

            //- Get width of output field
            virtual int width() const = 0;

            //- Set width of output field (and return old width)
            virtual int width(const int w) = 0;

            //- Get precision of output field
            virtual int precision() const = 0;

            //- Set precision of output field (and return old precision)
            virtual int precision(const int p) = 0;


    // Member Operators

        //- Return a non-const reference to const Ostream
        //  Needed for write functions where the stream argument is temporary:
        //  e.g. thing thisThing(OFstream("thingFileName")());
        Ostream& operator()() const
        {
            return const_cast<Ostream&>(*this);
        }


    // Housekeeping

        //- Access to indent level
        unsigned short& indentLevel() noexcept
        {
            return indentLevel_;
        }

        //- Access to indent size
        unsigned short& indentSize() noexcept
        {
            return indentSize_;
        }
};


// --------------------------------------------------------------------
// ------ Manipulators (not taking arguments)
// --------------------------------------------------------------------

//- An Ostream manipulator
typedef Ostream& (*OstreamManip)(Ostream&);

//- operator<< handling for manipulators without arguments
inline Ostream& operator<<(Ostream& os, OstreamManip f)
{
    return f(os);
}

//- operator<< handling for manipulators without arguments
inline Ostream& operator<<(Ostream& os, IOstreamManip f)
{
    f(os);
    return os;
}


//- Indent stream
inline Ostream& indent(Ostream& os)
{
    os.indent();
    return os;
}

//- Increment the indent level
inline Ostream& incrIndent(Ostream& os)
{
    os.incrIndent();
    return os;
}

//- Decrement the indent level
inline Ostream& decrIndent(Ostream& os)
{
    os.decrIndent();
    return os;
}


//- Flush stream
inline Ostream& flush(Ostream& os)
{
    os.flush();
    return os;
}


//- Add newline and flush stream
inline Ostream& endl(Ostream& os)
{
    os.endl();
    return os;
}


//- Write begin block group without a name
//  Increments indentation, adds newline.
inline Ostream& beginBlock(Ostream& os)
{
    os.beginBlock();
    return os;
}


//- Write end block group
//  Decrements indentation, adds newline.
inline Ostream& endBlock(Ostream& os)
{
    os.endBlock();
    return os;
}


//- Write end entry (';') followed by newline.
inline Ostream& endEntry(Ostream& os)
{
    os.endEntry();
    return os;
}


// Useful aliases for tab and newline characters
constexpr char tab = '\t';
constexpr char nl = '\n';


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
