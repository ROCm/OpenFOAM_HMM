/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2011 OpenCFD Ltd.
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

#include "conformalVoronoiMesh.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

bool Foam::conformalVoronoiMesh::insertSpecialisedFeaturePoint
(
    const featureEdgeMesh& feMesh,
    label ptI
)
{
    const labelList& pEds(feMesh.pointEdges()[ptI]);

    if (pEds.size() != 3)
    {
        // Only three edge specialisations available

        return false;
    }

    label nExternal = 0;
    label nInternal = 0;
    label nFlat = 0;
    label nOpen = 0;
    label nMultiple = 0;

    List<featureEdgeMesh::edgeStatus> allEdStat(pEds.size());

    forAll(pEds, i)
    {
        label edgeI = pEds[i];

        featureEdgeMesh::edgeStatus& eS = allEdStat[i];

        eS = feMesh.getEdgeStatus(edgeI);

        switch (eS)
        {
            case featureEdgeMesh::EXTERNAL:
            {
                nExternal++;
                break;
            }
            case featureEdgeMesh::INTERNAL:
            {
                nInternal++;
                break;
            }
            case featureEdgeMesh::FLAT:
            {
                nFlat++;
                break;
            }
            case featureEdgeMesh::OPEN:
            {
                nOpen++;
                break;
            }
            case featureEdgeMesh::MULTIPLE:
            {
                nMultiple++;
                break;
            }
            case featureEdgeMesh::NONE:
            {
                break;
            }
        }
    }

    if (nExternal == 2 && nInternal == 1)
    {
        Info<< "nExternal == 2 && nInternal == 1" << endl;

        return false;
    }
    else if (nExternal == 1 && nInternal == 2)
    {
        Info<< "nExternal == 1 && nInternal == 2" << endl;

        return false;
    }

    return false;
}


// ************************************************************************* //
