/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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

#include "lagrangianReconstructor.H"
#include "labelIOList.H"
#include "passiveParticleCloud.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lagrangianReconstructor::lagrangianReconstructor
(
    const fvMesh& mesh,
    const PtrList<fvMesh>& procMeshes,
    const PtrList<labelIOList>& faceProcAddressing,
    const PtrList<labelIOList>& cellProcAddressing
)
:
    mesh_(mesh),
    procMeshes_(procMeshes),
    faceProcAddressing_(faceProcAddressing),
    cellProcAddressing_(cellProcAddressing)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::lagrangianReconstructor::reconstructPositions
(
    const word& cloudName
) const
{
    passiveParticleCloud lagrangianPositions
    (
        mesh_,
        cloudName,
        IDLList<passiveParticle>()
    );

    forAll(procMeshes_, meshi)
    {
        const labelList& cellMap = cellProcAddressing_[meshi];
        const labelList& faceMap = faceProcAddressing_[meshi];

        Cloud<passiveParticle> lpi(procMeshes_[meshi], cloudName, false);

        forAllConstIters(lpi, iter)
        {
            const passiveParticle& ppi = *iter;

            const label mappedCell = cellMap[ppi.cell()];

            // Inverting sign if necessary and subtracting 1 from
            // faceProcAddressing
            const label mappedTetFace = mag(faceMap[ppi.tetFace()]) - 1;

            lagrangianPositions.append
            (
                new passiveParticle
                (
                    mesh_,
                    ppi.coordinates(),
                    mappedCell,
                    mappedTetFace,
                    ppi.procTetPt(mesh_, mappedCell, mappedTetFace)
                )
            );
        }
    }

    IOPosition<Cloud<passiveParticle>>(lagrangianPositions).write();

    // Force writing of "positions" too, if specified via the InfoSwitch
    if (particle::writeLagrangianPositions)
    {
        IOPosition<Cloud<passiveParticle>>
        (
            lagrangianPositions,
            cloud::geometryType::POSITIONS
        ).write();
    }

    return lagrangianPositions.size();
}


// ************************************************************************* //
