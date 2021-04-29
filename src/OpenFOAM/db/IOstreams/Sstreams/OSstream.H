/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
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
    Foam::OSstream

Description
    Generic output stream using a standard (STL) stream.

SourceFiles
    OSstreamI.H
    OSstream.C

\*---------------------------------------------------------------------------*/

#ifndef OSstream_H
#define OSstream_H

#include "Ostream.H"
#include "fileName.H"
#include <iostream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class OSstream Declaration
\*---------------------------------------------------------------------------*/

class OSstream
:
    public Ostream
{
    // Private Data

        fileName name_;

        std::ostream& os_;


public:

    // Generated Methods

        //- Copy construct
        OSstream(const OSstream&) = default;

        //- No copy assignment
        void operator=(const OSstream&) = delete;


    // Constructors

        //- Construct wrapper around std::ostream, set stream status
        //  Default stream options (ASCII, uncompressed)
        inline OSstream
        (
            std::ostream& os,
            const string& streamName,
            IOstreamOption streamOpt = IOstreamOption()
        );

        //- Construct wrapper around std::ostream, set stream status
        OSstream
        (
            std::ostream& os,
            const string& streamName,
            IOstreamOption::streamFormat fmt,
            IOstreamOption::versionNumber ver = IOstreamOption::currentVersion,
            IOstreamOption::compressionType cmp = IOstreamOption::UNCOMPRESSED
        )
        :
            OSstream(os, streamName, IOstreamOption(fmt, ver, cmp))
        {}


    // Member Functions

    // Characteristics

        //- Get the name of the stream.
        //  Useful for Fstream to remember the filename
        virtual const fileName& name() const
        {
            return name_;
        }

        //- Return stream name for modification
        virtual fileName& name()
        {
            return name_;
        }

        //- Get stream flags
        virtual ios_base::fmtflags flags() const;


    // Write Functions

            //- Write token to stream or otherwise handle it.
            //  \return false if the token type was not handled by this method
            virtual bool write(const token& tok);

            //- Write character
            virtual Ostream& write(const char c);

            //- Write character string
            virtual Ostream& write(const char* str);

            //- Write word
            virtual Ostream& write(const word& str);

            //- Write string (quoted)
            //  In the rare case that the string contains a final trailing
            //  backslash, it will be dropped to the appearance of an escaped
            //  double-quote.
            virtual Ostream& write(const string& str);

            //- Write std::string surrounded by quotes.
            //  Optional write without quotes.
            virtual Ostream& writeQuoted
            (
                const std::string& str,
                const bool quoted=true
            );

            //- Write int32_t
            virtual Ostream& write(const int32_t val);

            //- Write int64_t
            virtual Ostream& write(const int64_t val);

            //- Write floatScalar
            virtual Ostream& write(const floatScalar val);

            //- Write doubleScalar
            virtual Ostream& write(const doubleScalar val);

            //- Write binary block
            virtual Ostream& write(const char* data, std::streamsize count);

            //- Low-level raw binary output
            virtual Ostream& writeRaw
            (
                const char* data,
                std::streamsize count
            );

            //- Begin marker for low-level raw binary output.
            //  The count indicates the number of bytes for subsequent
            //  writeRaw calls.
            virtual bool beginRawWrite(std::streamsize count);

            //- End marker for low-level raw binary output.
            virtual bool endRawWrite();

            //- Add indentation characters
            virtual void indent();


        // Stream state functions

            //- Set stream flags
            virtual ios_base::fmtflags flags(const ios_base::fmtflags f);

            //- Flush stream
            virtual void flush();

            //- Add newline and flush stream
            virtual void endl();

            //- Get the current padding character
            virtual char fill() const;

            //- Set padding character for formatted field up to field width
            //  \return previous padding character
            virtual char fill(const char fillch);

            //- Get width of output field
            virtual int width() const;

            //- Set width of output field
            //  \return previous width
            virtual int width(const int w);

            //- Get precision of output field
            virtual int precision() const;

            //- Set precision of output field
            //  \return old precision
            virtual int precision(const int p);


        // STL stream

            //- Access to underlying std::ostream
            virtual std::ostream& stdStream()
            {
                return os_;
            }

            //- Const access to underlying std::ostream
            virtual const std::ostream& stdStream() const
            {
                return os_;
            }


        // Print

            //- Print stream description to Ostream
            virtual void print(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "OSstreamI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
