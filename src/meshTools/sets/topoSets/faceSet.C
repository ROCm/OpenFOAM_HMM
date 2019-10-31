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

#include "faceSet.H"
#include "polyMesh.H"
#include "mapPolyMesh.H"
#include "syncTools.H"
#include "mapDistributePolyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faceSet, 0);
    addToRunTimeSelectionTable(topoSet, faceSet, word);
    addToRunTimeSelectionTable(topoSet, faceSet, size);
    addToRunTimeSelectionTable(topoSet, faceSet, set);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faceSet::faceSet(const IOobject& obj)
:
    topoSet(obj, typeName)
{}


Foam::faceSet::faceSet
(
    const polyMesh& mesh,
    const word& name,
    readOption r,
    writeOption w
)
:
    topoSet(mesh, typeName, name, r, w)
{
    check(mesh.nFaces());
}


Foam::faceSet::faceSet
(
    const polyMesh& mesh,
    const word& name,
    const label size,
    writeOption w
)
:
    topoSet(mesh, name, size, w)
{}


Foam::faceSet::faceSet
(
    const polyMesh& mesh,
    const word& name,
    const topoSet& set,
    writeOption w
)
:
    topoSet(mesh, name, set, w)
{}


Foam::faceSet::faceSet
(
    const polyMesh& mesh,
    const word& name,
    const labelHashSet& labels,
    writeOption w
)
:
    topoSet(mesh, name, labels, w)
{}


Foam::faceSet::faceSet
(
    const polyMesh& mesh,
    const word& name,
    labelHashSet&& labels,
    writeOption w
)
:
    topoSet(mesh, name, std::move(labels), w)
{}


Foam::faceSet::faceSet
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

void Foam::faceSet::sync(const polyMesh& mesh)
{
    labelHashSet& labels = *this;

    // Convert to boolList
    // TBD: could change to using bitSet for the synchronization

    const label len = mesh.nFaces();

    boolList contents(len, false);

    for (const label facei : labels)
    {
        contents.set(facei);
    }

    syncTools::syncFaceList(mesh, contents, orEqOp<bool>());


    // Update labelHashSet

    for (label i=0; i < len; ++i)
    {
        if (contents.test(i))
        {
            labels.set(i);
        }
    }
}


Foam::label Foam::faceSet::maxSize(const polyMesh& mesh) const
{
    return mesh.nFaces();
}


void Foam::faceSet::updateMesh(const mapPolyMesh& morphMap)
{
    updateLabels(morphMap.reverseFaceMap());
}


void Foam::faceSet::distribute(const mapDistributePolyMesh& map)
{
    labelHashSet& labels = *this;

    boolList contents(map.nOldFaces(), false);

    for (const label facei : labels)
    {
        contents.set(facei);
    }

    map.distributeFaceData(contents);

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


void Foam::faceSet::writeDebug
(
    Ostream& os,
    const primitiveMesh& mesh,
    const label maxLen
) const
{
    topoSet::writeDebug(os, mesh.faceCentres(), maxLen);
}


// ************************************************************************* //
