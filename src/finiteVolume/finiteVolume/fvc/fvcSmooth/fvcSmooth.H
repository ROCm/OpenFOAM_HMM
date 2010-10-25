/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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

function
    Foam::fvc::smooth

Description
    Function that uses smoothData to apply spatial smoothing of a
    volume field using the FaceCellWave algorithm

SourceFiles
    fvcSmooth.C

\*---------------------------------------------------------------------------*/

#ifndef fvcSmooth_H
#define fvcSmooth_H

#include "volFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fvc
{
    void smooth
    (
        volScalarField& field,
        const scalar coeff
    );

    void spread
    (
        volScalarField& field,
        const volScalarField& alpha,
        const label nLayers
    );

    void sweep
    (
        volScalarField& field,
        const volScalarField& alpha,
        const label nLayers
    );
}
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
