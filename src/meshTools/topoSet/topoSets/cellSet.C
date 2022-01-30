/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2016-2018 OpenCFD Ltd.
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

#include "cellSet.H"
#include "mapPolyMesh.H"
#include "polyMesh.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"
#include "mapDistributePolyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(cellSet, 0);

addToRunTimeSelectionTable(topoSet, cellSet, word);
addToRunTimeSelectionTable(topoSet, cellSet, size);
addToRunTimeSelectionTable(topoSet, cellSet, set);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellSet::cellSet(const IOobject& obj)
:
    topoSet(obj, typeName)
{}


Foam::cellSet::cellSet
(
    const polyMesh& mesh,
    const word& name,
    readOption r,
    writeOption w
)
:
    topoSet(mesh, typeName, name, r, w)
{
    // Make sure set within valid range
    check(mesh.nCells());
}


Foam::cellSet::cellSet
(
    const polyMesh& mesh,
    const word& name,
    const label size,
    writeOption w
)
:
    topoSet(mesh, name, size, w)
{}


Foam::cellSet::cellSet
(
    const polyMesh& mesh,
    const word& name,
    const topoSet& set,
    writeOption w
)
:
    topoSet(mesh, name, set, w)
{}


Foam::cellSet::cellSet
(
    const polyMesh& mesh,
    const word& name,
    const labelHashSet& labels,
    writeOption w
)
:
    topoSet(mesh, name, labels, w)
{}


Foam::cellSet::cellSet
(
    const polyMesh& mesh,
    const word& name,
    labelHashSet&& labels,
    writeOption w
)
:
    topoSet(mesh, name, std::move(labels), w)
{}


Foam::cellSet::cellSet
(
    const polyMesh& mesh,
    const word& name,
    const labelUList& labels,
    writeOption w
)
:
    topoSet(mesh, name, labels, w)
{}


// Database constructors (for when no mesh available)
Foam::cellSet::cellSet
(
    const Time& runTime,
    const word& name,
    readOption r,
    writeOption w
)
:
    topoSet
    (
        findIOobject(runTime, name, r, w),
        typeName
    )
{}


Foam::cellSet::cellSet
(
    const Time& runTime,
    const word& name,
    const label size,
    writeOption w
)
:
    topoSet
    (
        findIOobject(runTime, name, IOobject::NO_READ, w),
        size
    )
{}


Foam::cellSet::cellSet
(
    const Time& runTime,
    const word& name,
    const labelHashSet& labels,
    writeOption w
)
:
    topoSet
    (
        findIOobject(runTime, name, IOobject::NO_READ, w),
        labels
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::cellSet::maxSize(const polyMesh& mesh) const
{
    return mesh.nCells();
}


void Foam::cellSet::updateMesh(const mapPolyMesh& morphMap)
{
    updateLabels(morphMap.reverseCellMap());
}


void Foam::cellSet::distribute(const mapDistributePolyMesh& map)
{
    labelHashSet& labels = *this;

    boolList contents(map.nOldCells(), false);

    for (const label celli : labels)
    {
        contents.set(celli);
    }

    map.distributeCellData(contents);

    // The new length
    const label len = contents.size();

    // Count
    label n = 0;
    for (label i=0; i < len; ++i)
    {
        if (contents.test(i))
        {
            ++n;
        }
    }

    // Update labelHashSet

    labels.clear();
    labels.resize(2*n);

    for (label i=0; i < len; ++i)
    {
        if (contents.test(i))
        {
            labels.set(i);
        }
    }
}


void Foam::cellSet::writeDebug
(
    Ostream& os,
    const primitiveMesh& mesh,
    const label maxLen
) const
{
    topoSet::writeDebug(os, mesh.cellCentres(), maxLen);
}


// ************************************************************************* //
