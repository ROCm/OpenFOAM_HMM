/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
     \\/     M anipulation  |
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
    Foam::vtk::legacyRawFormatter

Description
    Binary output for the VTK legacy format, always written as big-endian
    and with 32-bit integers.

    This format should never be used for OpenFOAM with 64-bit label sizes.

SourceFiles
    foamVtkLegacyRawFormatter.C

\*---------------------------------------------------------------------------*/

#ifndef foamVtkLegacyRawFormatter_H
#define foamVtkLegacyRawFormatter_H

#include "foamVtkFormatter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace vtk
{

/*---------------------------------------------------------------------------*\
                       Class legacyRawFormatter Declaration
\*---------------------------------------------------------------------------*/

class legacyRawFormatter
:
    public formatter
{
    // Private Data Members

        static const char* legacyName_;
        static const outputOptions opts_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        legacyRawFormatter(const legacyRawFormatter&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const legacyRawFormatter&) = delete;


protected:

    // Protected Member Functions

        //- Write
        void write(const char* s, std::streamsize n);


public:

    // Constructors

        //- Construct and attach to an output stream
        legacyRawFormatter(std::ostream& os);


    //- Destructor
    virtual ~legacyRawFormatter();


    // Member Functions

        //- The output is LEGACY_BINARY.
        virtual const outputOptions& opts() const;

        //- Name for the legacy binary output type ("BINARY")
        virtual const char* name() const;

        //- Name for the XML append encoding (unused)
        //  Currently identical to name(), but do not rely on this.
        virtual const char* encoding() const;


        //- Write leading size - a no-op for legacy binary output
        virtual void writeSize(const uint64_t ignored);

        virtual void write(const uint8_t val);
        virtual void write(const label val);
        virtual void write(const float val);
        virtual void write(const double val);

        //- Write a newline to the output
        virtual void flush();

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
