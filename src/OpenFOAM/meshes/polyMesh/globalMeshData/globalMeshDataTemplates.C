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

#include "globalMeshData.H"
#include "polyMesh.H"
#include "mapDistribute.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class CombineOp>
void Foam::globalMeshData::syncPointData
(
    List<Type>& pointData,
    const labelListList& slaves,
    const mapDistribute& slavesMap,
    const CombineOp& cop
) const
{
    if (pointData.size() != mesh_.nPoints())
    {
        FatalErrorIn("globalMeshData::syncPointData(..)")
            << "Number of elements in data:" << pointData.size()
            << " differs from number of points in mesh:" << mesh_.nPoints()
            << abort(FatalError);
    }

    const indirectPrimitivePatch& cpp = coupledPatch();
    const labelList& meshPoints = cpp.meshPoints();

    // Copy mesh (point)data to coupled patch (point)data
    Field<Type> cppFld(slavesMap.constructSize());
    forAll(meshPoints, patchPointI)
    {
        cppFld[patchPointI] = pointData[meshPoints[patchPointI]];
    }

    // Pull slave data onto master
    slavesMap.distribute(cppFld);

    // Combine master data with slave data
    forAll(slaves, patchPointI)
    {
        const labelList& slavePoints = slaves[patchPointI];

        // Combine master with slave data
        forAll(slavePoints, i)
        {
            cop(cppFld[patchPointI], cppFld[slavePoints[i]]);
        }
        // Copy result back to slave slots
        forAll(slavePoints, i)
        {
            cppFld[slavePoints[i]] = cppFld[patchPointI];
        }
    }

    // Push master data back to slaves
    slavesMap.reverseDistribute(meshPoints.size(), cppFld);

    // Update mesh (point)data from coupled patch (point)data
    forAll(meshPoints, patchPointI)
    {
        pointData[meshPoints[patchPointI]] = cppFld[patchPointI];
    }
}


template<class Type, class CombineOp>
void Foam::globalMeshData::syncPointData
(
    List<Type>& pointData,
    const CombineOp& cop
) const
{
    const labelListList& slaves = globalPointSlaves();
    const mapDistribute& map = globalPointSlavesMap();

    syncPointData
    (
        pointData,
        slaves,
        map,
        cop
    );
}


template<class Type, class CombineOp>
void Foam::globalMeshData::syncPointAllData
(
    List<Type>& pointData,
    const CombineOp& cop
) const
{
    syncPointData
    (
        pointData,
        globalPointAllSlaves(),
        globalPointAllSlavesMap(),
        cop
    );
}


// ************************************************************************* //
