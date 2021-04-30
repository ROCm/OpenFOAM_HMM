/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::faSolution

Description
    Selector class for finite area solution.
    faMesh is derived from faSolution so that all fields have access to the
    faSolution from the mesh reference they hold.

Author
    Zeljko Tukovic, FMENA
    Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#ifndef faSolution_H
#define faSolution_H

#include "solution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class faSolution Declaration
\*---------------------------------------------------------------------------*/

class faSolution
:
    public solution
{
public:

    //- No copy construct
    faSolution(const faSolution&) = delete;

    //- No copy assignment
    void operator=(const faSolution&) = delete;


    // Constructors

        //- Construct for objectRegistry, dictionary name and optional
        //- fallback dictionary content (for a NO_READ or missing file)
        //  A null dictionary pointer is treated like an empty dictionary.
        faSolution
        (
            const objectRegistry& obr,
            const word& dictName,
            const dictionary* fallback = nullptr
        )
        :
            solution(obr, dictName, fallback)
        {}

        //- Construct for objectRegistry, and optional
        //- fallback dictionary content (for a NO_READ or missing file)
        //  A null dictionary pointer is treated like an empty dictionary.
        explicit faSolution
        (
            const objectRegistry& obr,
            const dictionary* fallback = nullptr
        )
        :
            faSolution(obr, "faSolution", fallback)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
