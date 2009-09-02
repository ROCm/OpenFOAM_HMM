/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "mathConstants.H"
#include "universalConstants.H"
#include "electromagneticConstants.H"
#include "physicoChemicalConstants.H"

#include "dimensionedConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const char* Foam::constant::physicoChemical::group = "physicoChemical";


const Foam::dimensionedScalar Foam::constant::physicoChemical::R
(
    dimensionedConstant
    (
        group,
        "R",
        dimensionedScalar
        (
            "R",
            NA*k
        )
    )
);


const Foam::dimensionedScalar Foam::constant::physicoChemical::F
(
    dimensionedConstant
    (
        group,
        "F",
        dimensionedScalar
        (
            "F",
            NA*constant::electromagnetic::e
        )
    )
);


const Foam::dimensionedScalar Foam::constant::physicoChemical::sigma
(
    dimensionedConstant
    (
        group,
        "sigma",
        dimensionedScalar
        (
            "sigma",
            dimensionedScalar("C", dimless, sqr(constant::math::pi)/60.0)
           *pow4(k)/(pow3(constant::universal::hr)*sqr(constant::universal::c))
        )
    )
);


const Foam::dimensionedScalar Foam::constant::physicoChemical::b
(
    dimensionedConstant
    (
        group,
        "b",
        dimensionedScalar
        (
            "b",
            (constant::universal::h*constant::universal::c/k)
           /dimensionedScalar("C", dimless, 4.965114231)
        )
    )
);


const Foam::dimensionedScalar Foam::constant::physicoChemical::c1
(
    dimensionedConstant
    (
        group,
        "c1",
        dimensionedScalar
        (
            "c1",
            dimensionedScalar("C", dimless, constant::math::twoPi)
           *constant::universal::h*sqr(constant::universal::c)
        )
    )
);


const Foam::dimensionedScalar Foam::constant::physicoChemical::c2
(
    dimensionedConstant
    (
        group,
        "c2",
        dimensionedScalar
        (
            "c2",
            constant::universal::h*constant::universal::c/k
        )
    )
);


// ************************************************************************* //
