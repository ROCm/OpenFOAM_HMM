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

#include "rampHoldFall.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(rampHoldFall, 0);
addToRunTimeSelectionTable(relaxationModel, rampHoldFall, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

rampHoldFall::rampHoldFall
(
    const dictionary& relaxationDict,
    const conformalVoronoiMesh& cvMesh
)
:
    relaxationModel(typeName, relaxationDict, cvMesh),
    rampStartRelaxation_(readScalar(coeffDict().lookup("rampStartRelaxation"))),
    holdRelaxation_(readScalar(coeffDict().lookup("holdRelaxation"))),
    fallEndRelaxation_(readScalar(coeffDict().lookup("fallEndRelaxation"))),
    rampEndFraction_(readScalar(coeffDict().lookup("rampEndFraction"))),
    fallStartFraction_(readScalar(coeffDict().lookup("fallStartFraction"))),
    rampGradient_((holdRelaxation_ - rampStartRelaxation_)/(rampEndFraction_)),
    fallGradient_
    (
        (fallEndRelaxation_ - holdRelaxation_)/(1 - fallStartFraction_)
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar rampHoldFall::relaxation()
{
    scalar t = cvMesh_.time().timeOutputValue();

    scalar tStart = cvMesh_.time().startTime().value();
    scalar tEnd = cvMesh_.time().endTime().value();
    scalar tSpan = tEnd - tStart;

    if (tSpan < VSMALL)
    {
        return rampStartRelaxation_;
    }

    if (t - tStart < rampEndFraction_*tSpan)
    {
        // Ramp

        return rampGradient_*((t - tStart)/tSpan) + rampStartRelaxation_;
    }
    else if (t - tStart > fallStartFraction_*tSpan)
    {
        // Fall

        return
            fallGradient_*((t - tStart)/tSpan)
          + fallEndRelaxation_ - fallGradient_;
    }
    else
    {
        //Hold

        return holdRelaxation_;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
