/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "piecewiseLinearRamp.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(piecewiseLinearRamp, 0);
addToRunTimeSelectionTable(faceAreaWeightModel, piecewiseLinearRamp, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

piecewiseLinearRamp::piecewiseLinearRamp
(
    const dictionary& faceAreaWeightDict,
    const conformalVoronoiMesh& cvMesh
)
:
    faceAreaWeightModel(typeName, faceAreaWeightDict, cvMesh),
    lAF_(readScalar(coeffDict().lookup("lowerAreaFraction"))),
    uAF_(readScalar(coeffDict().lookup("upperAreaFraction")))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar piecewiseLinearRamp::faceAreaWeight(scalar faceAreaFraction) const
{
    if (faceAreaFraction < lAF_)
    {
        return 0;
    }
    else if (faceAreaFraction < uAF_)
    {
        return faceAreaFraction/((uAF_ - lAF_)) - 1/((uAF_/lAF_) - 1);
    }
    else
    {
        return 1;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
