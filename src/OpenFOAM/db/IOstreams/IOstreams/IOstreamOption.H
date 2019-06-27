/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::IOstreamOption

Description
    The IOstreamOption is a simple container for options an IOstream
    can normally have. For example, ascii/binary, uncompressed/compressed, ...

SourceFiles
    IOstreamOption.C

\*---------------------------------------------------------------------------*/

#ifndef IOstreamOption_H
#define IOstreamOption_H

#include "scalar.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class Ostream;
template<class EnumType> class Enum;

/*---------------------------------------------------------------------------*\
                       Class IOstreamOption Declaration
\*---------------------------------------------------------------------------*/

class IOstreamOption
{
public:

    // Public Data Types

        //- Data format (ascii | binary)
        enum streamFormat : char
        {
            ASCII,              //!< "ascii"
            BINARY              //!< "binary"
        };

        //- Compression treatment (UNCOMPRESSED | COMPRESSED)
        enum compressionType : char
        {
            UNCOMPRESSED = 0,   //!< compression = false
            COMPRESSED          //!< compression = true
        };


        //- Representation of a major/minor version number
        class versionNumber
        {
            //- The combined major/version number.
            short number_;

        public:

        // Constructors

            //- Construct from major, number
            constexpr versionNumber(int major, int minor) noexcept
            :
                number_(10*major + (minor % 10))
            {}

            //- Construct from floating-point version number
            explicit constexpr versionNumber(const float ver) noexcept
            :
                number_(10*ver + 0.001) // Allow some rounding
            {}

            //- Construct from Istream by reading in a float.
            //  Non-explicit for convenience
            versionNumber(Istream& is)
            :
                versionNumber(readFloat(is))
            {}

            //- Construct by parsing string "major.minor"
            explicit versionNumber(const std::string& verNum)
            :
                versionNumber(readFloat(verNum))
            {}


        // Member Functions

            //- The canonical major/minor pair as an integer value.
            inline int canonical() noexcept
            {
                return number_;
            }

            //- Return the major version number.
            inline int getMajor() const noexcept
            {
                return int(number_ / 10);
            }

            //- Return the minor version number
            inline int getMinor() const noexcept
            {
                return int(number_ % 10);
            }

            //- A string representation of major.minor
            std::string str() const
            {
                return
                    std::to_string(getMajor())
                  + '.'
                  + std::to_string(getMinor());
            }


        // Member Operators

            //- Version number equality
            bool operator==(const versionNumber& rhs) const noexcept
            {
                return number_ == rhs.number_;
            }

            //- Version number inequality
            bool operator!=(const versionNumber& rhs) const noexcept
            {
                return number_ != rhs.number_;
            }

            //- Version number older than rhs
            bool operator<(const versionNumber& rhs) const noexcept
            {
                return number_ < rhs.number_;
            }

            //- Version number is the same or older than rhs
            bool operator<=(const versionNumber& rhs) const noexcept
            {
                return number_ <= rhs.number_;
            }

            //- Version number newer than rhs
            bool operator>(const versionNumber& rhs) const noexcept
            {
                return number_ > rhs.number_;
            }

            //- Version number same or newer than rhs
            bool operator>=(const versionNumber& rhs) const noexcept
            {
                return number_ >= rhs.number_;
            }
        };


    // Public Static Data

        //- Stream format names (ascii, binary)
        static const Enum<streamFormat> formatNames;

        //- The original version number
        static const versionNumber originalVersion;

        //- The current version number
        static const versionNumber currentVersion;


private:

    // Private Data

    // NB: ordered with adjacent enums to minimize gaps

        //- Stream version number (eg, 2.0 for current dictionary format)
        versionNumber version_;

        //- Format: (ascii | binary)
        streamFormat format_;

        //- Compression: (on | off)
        compressionType compression_;


public:

    // Constructors

        //- Construct null. (default: ASCII, uncompressed, currentVersion)
        IOstreamOption() noexcept
        :
            version_(currentVersion),
            format_(ASCII),
            compression_(compressionType::UNCOMPRESSED)
        {}

        //- Construct with format. (default: uncompressed, currentVersion)
        explicit IOstreamOption(streamFormat format) noexcept
        :
            version_(currentVersion),
            format_(format),
            compression_(compressionType::UNCOMPRESSED)
        {}

        //- Construct with format and compression, optionally with version.
        IOstreamOption
        (
            streamFormat format,
            compressionType compression,
            versionNumber version=currentVersion
        ) noexcept
        :
            version_(version),
            format_(format),
            compression_(compression)
        {}

        //- Construct with format, version, compression
        IOstreamOption
        (
            streamFormat format,
            versionNumber version,
            compressionType compression
        ) noexcept
        :
            version_(version),
            format_(format),
            compression_(compression)
        {}


    // Static Member Functions

        //- The stream format enum corresponding to the string
        //  Expected "ascii", "binary"
        static streamFormat formatEnum(const word& formatName);

        //- The compression enum corresponding to the string
        //  Expected "true", "false", "on", "off", etc.
        static compressionType compressionEnum(const word& compName);


    // Member Functions

        //- Get the current stream format
        streamFormat format() const noexcept
        {
            return format_;
        }

        //- Set the stream format
        //  \return the previous value
        streamFormat format(const streamFormat format) noexcept
        {
            streamFormat old(format_);
            format_ = format;
            return old;
        }

        //- Set the stream format, from string value
        //  \return the previous value
        streamFormat format(const word& formatName)
        {
            streamFormat old(format_);
            format_ = formatEnum(formatName);
            return old;
        }

        //- Get the stream compression
        compressionType compression() const noexcept
        {
            return compression_;
        }

        //- Set the stream compression
        //  \return the previous value
        compressionType compression(const compressionType comp) noexcept
        {
            compressionType old(compression_);
            compression_ = comp;
            return old;
        }

        //- Set the stream compression, from string value.
        //  \return the previous value
        compressionType compression(const word& compressionName)
        {
            compressionType old(compression_);
            compression_ = compressionEnum(compressionName);
            return old;
        }

        //- Get the stream version
        versionNumber version() const noexcept
        {
            return version_;
        }

        //- Set the stream version
        //  \return the previous value
        versionNumber version(const versionNumber verNum) noexcept
        {
            versionNumber old(version_);
            version_ = verNum;
            return old;
        }
};


//- Output the format as text string (ascii | binary)
Ostream& operator<<(Ostream& os, const IOstreamOption::streamFormat& sf);

//- Output the version as major.minor
Ostream& operator<<(Ostream& os, const IOstreamOption::versionNumber& vn);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
