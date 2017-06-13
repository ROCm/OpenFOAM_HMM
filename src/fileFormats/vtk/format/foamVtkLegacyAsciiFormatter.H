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
    Foam::vtk::legacyAsciiFormatter

Description
    Formatting as per Foam::vtk::asciiFormatter, but with
    naming for legacy output.

SourceFiles
    foamVtkLegacyAsciiFormatter.C

\*---------------------------------------------------------------------------*/

#ifndef foamVtkLegacyAsciiFormatter_H
#define foamVtkLegacyAsciiFormatter_H

#include "foamVtkAsciiFormatter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace vtk
{

/*---------------------------------------------------------------------------*\
                    Class legacyAsciiFormatter Declaration
\*---------------------------------------------------------------------------*/

class legacyAsciiFormatter
:
    public asciiFormatter
{
    // Private Data Members

        static const char* legacyName_;
        static const outputOptions opts_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        legacyAsciiFormatter(const legacyAsciiFormatter&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const legacyAsciiFormatter&) = delete;


public:

    // Constructors

        //- Construct and attach to an output stream, use default precision
        legacyAsciiFormatter(std::ostream& os);

        //- Construct and attach to an output stream, use specified precision
        legacyAsciiFormatter(std::ostream& os, unsigned precision);


    //- Destructor
    virtual ~legacyAsciiFormatter();


    // Member Functions

        //- The output is LEGACY_ASCII.
        virtual const outputOptions& opts() const;

        //- Name for the legacy ascii output type ("ASCII")
        virtual const char* name() const;

        //- Name for the XML append encoding (unused).
        //  Currently identical to name(), but do not rely on this.
        virtual const char* encoding() const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
