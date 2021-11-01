/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
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
    Foam::prefixOSstream

Description
    Version of OSstream that prints a prefix on each line.

    This is useful for running in parallel as it allows the processor number
    to be automatically prepended to each message line.

SourceFiles
    prefixOSstream.C

\*---------------------------------------------------------------------------*/

#ifndef prefixOSstream_H
#define prefixOSstream_H

#include "OSstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class prefixOSstream Declaration
\*---------------------------------------------------------------------------*/

class prefixOSstream
:
    public OSstream
{
    // Private Data

        bool printPrefix_;

        string prefix_;


    // Private Member Functions

        inline void checkWritePrefix();


public:

    // Constructors

        //- Construct and set stream status
        prefixOSstream
        (
            std::ostream& os,
            const string& streamName,
            IOstreamOption streamOpt = IOstreamOption()
        );

        //- Construct and set stream status
        prefixOSstream
        (
            std::ostream& os,
            const string& streamName,
            IOstreamOption::streamFormat fmt,
            IOstreamOption::versionNumber ver = IOstreamOption::currentVersion,
            IOstreamOption::compressionType cmp = IOstreamOption::UNCOMPRESSED
        )
        :
            prefixOSstream(os, streamName, IOstreamOption(fmt, ver, cmp))
        {}


    // Member Functions

        // Enquiry

            //- Return the stream prefix
            const string& prefix() const noexcept
            {
                return prefix_;
            }

            //- Return non-const access to the stream prefix
            string& prefix() noexcept
            {
                return prefix_;
            }


        // Write Functions

            //- Write token to stream or otherwise handle it.
            //  \return false if the token type was not handled by this method
            virtual bool write(const token& tok);

            //- Write character
            virtual Ostream& write(const char c);

            //- Write character string
            virtual Ostream& write(const char* str);

            //- Write word
            virtual Ostream& write(const word& val);

            //- Write string
            virtual Ostream& write(const string& val);

            //- Write std::string surrounded by quotes.
            //  Optional write without quotes.
            virtual Ostream& writeQuoted
            (
                const std::string& val,
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
            virtual Ostream& write(const char* buf, std::streamsize count);

            //- Add indentation characters
            virtual void indent();


    // Print

        //- Print stream description to Ostream
        virtual void print(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
