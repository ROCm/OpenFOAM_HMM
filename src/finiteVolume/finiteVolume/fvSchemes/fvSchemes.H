/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::fvSchemes

Description
    Selector class for finite volume differencing schemes.
    fvMesh is derived from fvSchemes so that all fields have access to the
    fvSchemes from the mesh reference they hold.

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef fvSchemes_H
#define fvSchemes_H

#include "schemesLookup.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fvSchemes Declaration
\*---------------------------------------------------------------------------*/

class fvSchemes
:
    public schemesLookup
{
public:

    // Constructors

        //- No copy construct
        fvSchemes(const fvSchemes&) = delete;

        //- No copy assignment
        void operator=(const fvSchemes&) = delete;

        //- Construct for objectRegistry, dictionary name and optional
        //- fallback dictionary content (for a NO_READ or missing file)
        //  A null dictionary pointer is treated like an empty dictionary.
        fvSchemes
        (
            const objectRegistry& obr,
            const word& dictName,
            const dictionary* fallback = nullptr
        )
        :
            schemesLookup(obr, dictName, fallback)
        {}

        //- Construct for objectRegistry, and optional
        //- fallback dictionary content (for a NO_READ or missing file)
        //  A null dictionary pointer is treated like an empty dictionary.
        explicit fvSchemes
        (
            const objectRegistry& obr,
            const dictionary* fallback = nullptr
        )
        :
            fvSchemes(obr, "fvSchemes", fallback)
        {}

        //- Construct for objectRegistry, and
        //- fallback dictionary content (for a NO_READ or missing file)
        fvSchemes(const objectRegistry& obr, const dictionary& dict)
        :
            fvSchemes(obr, "fvSchemes", &dict)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
