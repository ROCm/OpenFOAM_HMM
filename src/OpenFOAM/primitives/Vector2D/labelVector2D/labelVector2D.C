/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "labelVector2D.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const Foam::labelVector2D::vsType::typeName = "labelVector2D";

template<>
const char* const Foam::labelVector2D::vsType::componentNames[] = { "x", "y" };

template<>
const Foam::labelVector2D Foam::labelVector2D::vsType::zero
(
    labelVector2D::uniform(0)
);

template<>
const Foam::labelVector2D Foam::labelVector2D::vsType::one
(
    labelVector2D::uniform(1)
);

template<>
const Foam::labelVector2D Foam::labelVector2D::vsType::max
(
    labelVector2D::uniform(labelMax)
);

template<>
const Foam::labelVector2D Foam::labelVector2D::vsType::min
(
    labelVector2D::uniform(-labelMax)
);

template<>
const Foam::labelVector2D Foam::labelVector2D::vsType::rootMax
(
    labelVector2D::uniform(sqrt(scalar(labelMax)))
);

template<>
const Foam::labelVector2D Foam::labelVector2D::vsType::rootMin
(
    labelVector2D::uniform(-sqrt(scalar(labelMax)))
);


// ************************************************************************* //
