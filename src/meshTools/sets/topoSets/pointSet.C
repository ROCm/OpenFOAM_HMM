/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "pointSet.H"
#include "mapPolyMesh.H"
#include "polyMesh.H"
#include "syncTools.H"
#include "mapDistributePolyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pointSet, 0);
    addToRunTimeSelectionTable(topoSet, pointSet, word);
    addToRunTimeSelectionTable(topoSet, pointSet, size);
    addToRunTimeSelectionTable(topoSet, pointSet, set);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::pointSet::pointSet(const IOobject& obj)
:
    topoSet(obj, typeName)
{}


Foam::pointSet::pointSet
(
    const polyMesh& mesh,
    const word& name,
    readOption r,
    writeOption w
)
:
    topoSet(mesh, typeName, name, r, w)
{
    check(mesh.nPoints());
}


Foam::pointSet::pointSet
(
    const polyMesh& mesh,
    const word& name,
    const label size,
    writeOption w
)
:
    topoSet(mesh, name, size, w)
{}


Foam::pointSet::pointSet
(
    const polyMesh& mesh,
    const word& name,
    const topoSet& set,
    writeOption w
)
:
    topoSet(mesh, name, set, w)
{}


Foam::pointSet::pointSet
(
    const polyMesh& mesh,
    const word& name,
    const labelHashSet& labels,
    writeOption w
)
:
    topoSet(mesh, name, labels, w)
{}


Foam::pointSet::pointSet
(
    const polyMesh& mesh,
    const word& name,
    labelHashSet&& labels,
    writeOption w
)
:
    topoSet(mesh, name, std::move(labels), w)
{}


Foam::pointSet::pointSet
(
    const polyMesh& mesh,
    const word& name,
    const labelUList& labels,
    writeOption w
)
:
    topoSet(mesh, name, labels, w)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pointSet::sync(const polyMesh& mesh)
{
    labelHashSet& labels = *this;

    // Convert to boolList
    // TBD: could change to using bitSet for the synchronization

    const label len = mesh.nPoints();

    boolList contents(len, false);

    for (const label pointi : labels)
    {
        contents.set(pointi);
    }

    // The nullValue = 'false'
    syncTools::syncPointList(mesh, contents, orEqOp<bool>(), false);


    // Update labelHashSet

    for (label pointi=0; pointi < len; ++pointi)
    {
        if (contents.test(pointi))
        {
            labels.set(pointi);
        }
    }
}


Foam::label Foam::pointSet::maxSize(const polyMesh& mesh) const
{
    return mesh.nPoints();
}


void Foam::pointSet::updateMesh(const mapPolyMesh& morphMap)
{
    updateLabels(morphMap.reversePointMap());
}


void Foam::pointSet::distribute(const mapDistributePolyMesh& map)
{
    labelHashSet& labels = *this;

    boolList contents(map.nOldPoints(), false);

    for (const label pointi : labels)
    {
        contents.set(pointi);
    }

    map.distributePointData(contents);

    // The new length
    const label len = contents.size();

    // Count - as per BitOps::count(contents)
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


void Foam::pointSet::writeDebug
(
    Ostream& os,
    const primitiveMesh& mesh,
    const label maxLen
) const
{
    topoSet::writeDebug(os, mesh.points(), maxLen);
}


// ************************************************************************* //
