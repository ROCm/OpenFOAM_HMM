/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "uniformValue.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(uniformValue, 0);
    addToRunTimeSelectionTable
    (
        surfaceCellSizeFunction,
        uniformValue,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uniformValue::uniformValue
(
    const dictionary& cellSizeFunctionDict,
    const searchableSurface& surface
)
:
    surfaceCellSizeFunction(typeName, cellSizeFunctionDict, surface),
    surfaceCellSize_(readScalar(coeffsDict().lookup("surfaceCellSize")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::scalar& Foam::uniformValue::surfaceSize(const label index) const
{
    return surfaceCellSize_;
}


const Foam::scalar& Foam::uniformValue::refineSurfaceSize(const label index)
{
    surfaceCellSize_ *= refinementFactor_;

    return surfaceCellSize_;
}


Foam::scalar Foam::uniformValue::interpolate
(
    const point& pt,
    const label index
) const
{
    return surfaceCellSize_;
}


void Foam::uniformValue::recalculateInterpolation() const
{}


// ************************************************************************* //
