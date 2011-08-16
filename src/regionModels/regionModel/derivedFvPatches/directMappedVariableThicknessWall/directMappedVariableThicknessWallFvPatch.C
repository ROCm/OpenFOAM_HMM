/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "directMappedVariableThicknessWallFvPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "regionModel1D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(directMappedVariableThicknessWallFvPatch, 0);
    addToRunTimeSelectionTable
    (
        fvPatch,
        directMappedVariableThicknessWallFvPatch,
        polyPatch
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::directMappedVariableThicknessWallFvPatch::
makeDeltaCoeffs(scalarField& dc) const
{
    const directMappedVariableThicknessWallPolyPatch& pp =
        refCast<const directMappedVariableThicknessWallPolyPatch>
        (
            patch()
        );

    const directMappedPatchBase& mpp = refCast<const directMappedPatchBase>
    (
        patch()
    );

    const polyMesh& nbrMesh = mpp.sampleMesh();

    typedef regionModels::regionModel1D modelType;

    const modelType& region1D =
        nbrMesh.objectRegistry::lookupObject<modelType>
        (
            "thermoBaffleProperties"
        );

    dc = 2.0/(pp.thickness()/region1D.nLayers());
}


// ************************************************************************* //
