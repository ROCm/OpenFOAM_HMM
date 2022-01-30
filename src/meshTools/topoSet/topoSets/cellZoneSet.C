/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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

#include "cellZoneSet.H"
#include "mapPolyMesh.H"
#include "polyMesh.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cellZoneSet, 0);
    addToRunTimeSelectionTable(topoSet, cellZoneSet, word);
    addToRunTimeSelectionTable(topoSet, cellZoneSet, size);
    addToRunTimeSelectionTable(topoSet, cellZoneSet, set);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::cellZoneSet::updateSet()
{
    labelList order(sortedOrder(addressing_));
    inplaceReorder(order, addressing_);

    cellSet::clearStorage();
    cellSet::resize(2*addressing_.size());
    cellSet::set(addressing_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellZoneSet::cellZoneSet
(
    const polyMesh& mesh,
    const word& name,
    readOption r,
    writeOption w
)
:
    cellSet(mesh, name, 1024),  // do not read cellSet
    mesh_(mesh),
    addressing_()
{
    const cellZoneMesh& cellZones = mesh.cellZones();
    label zoneID = cellZones.findZoneID(name);

    if
    (
        (r == IOobject::MUST_READ)
     || (r == IOobject::MUST_READ_IF_MODIFIED)
     || (r == IOobject::READ_IF_PRESENT && zoneID != -1)
    )
    {
        const cellZone& fz = cellZones[zoneID];
        addressing_ = fz;
    }

    updateSet();

    check(mesh.nCells());
}


Foam::cellZoneSet::cellZoneSet
(
    const polyMesh& mesh,
    const word& name,
    const label size,
    writeOption w
)
:
    cellSet(mesh, name, size, w),
    mesh_(mesh),
    addressing_()
{
    updateSet();
}


Foam::cellZoneSet::cellZoneSet
(
    const polyMesh& mesh,
    const word& name,
    const topoSet& set,
    writeOption w
)
:
    cellSet(mesh, name, set.size(), w),
    mesh_(mesh),
    addressing_(refCast<const cellZoneSet>(set).addressing())
{
    updateSet();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cellZoneSet::invert(const label maxLen)
{
    // Count
    label n = 0;

    for (label celli = 0; celli < maxLen; ++celli)
    {
        if (!found(celli))
        {
            ++n;
        }
    }

    // Fill
    addressing_.setSize(n);
    n = 0;

    for (label celli = 0; celli < maxLen; ++celli)
    {
        if (!found(celli))
        {
            addressing_[n] = celli;
            ++n;
        }
    }

    updateSet();
}


void Foam::cellZoneSet::subset(const topoSet& set)
{
    DynamicList<label> newAddressing(addressing_.size());

    const cellZoneSet& zoneSet = refCast<const cellZoneSet>(set);

    for (const label celli : zoneSet.addressing())
    {
        if (found(celli))
        {
            newAddressing.append(celli);
        }
    }

    addressing_.transfer(newAddressing);
    updateSet();
}


void Foam::cellZoneSet::addSet(const topoSet& set)
{
    DynamicList<label> newAddressing(addressing_);

    const cellZoneSet& zoneSet = refCast<const cellZoneSet>(set);

    for (const label celli : zoneSet.addressing())
    {
        if (!found(celli))
        {
            newAddressing.append(celli);
        }
    }

    addressing_.transfer(newAddressing);
    updateSet();
}


void Foam::cellZoneSet::subtractSet(const topoSet& set)
{
    DynamicList<label> newAddressing(addressing_.size());

    const cellZoneSet& zoneSet = refCast<const cellZoneSet>(set);

    for (const label celli : addressing_)
    {
        if (!zoneSet.found(celli))
        {
            // Not found in zoneSet so add
            newAddressing.append(celli);
        }
    }

    addressing_.transfer(newAddressing);
    updateSet();
}


void Foam::cellZoneSet::sync(const polyMesh& mesh)
{
    cellSet::sync(mesh);

    // Take over contents of cellSet into addressing.
    addressing_ = sortedToc();
    updateSet();
}


Foam::label Foam::cellZoneSet::maxSize(const polyMesh& mesh) const
{
    return mesh.nCells();
}


bool Foam::cellZoneSet::writeObject
(
    IOstreamOption streamOpt,
    const bool valid
) const
{
    // Write shadow cellSet
    word oldTypeName = typeName;
    const_cast<word&>(type()) = cellSet::typeName;
    bool ok = cellSet::writeObject(streamOpt, valid);
    const_cast<word&>(type()) = oldTypeName;

    // Modify cellZone
    cellZoneMesh& cellZones = const_cast<polyMesh&>(mesh_).cellZones();
    label zoneID = cellZones.findZoneID(name());

    if (zoneID == -1)
    {
        zoneID = cellZones.size();

        cellZones.setSize(zoneID+1);
        cellZones.set
        (
            zoneID,
            new cellZone
            (
                name(),
                addressing_,
                zoneID,
                cellZones
            )
        );
    }
    else
    {
        cellZones[zoneID] = addressing_;
    }
    cellZones.clearAddressing();

    return ok && cellZones.write(valid);
}


void Foam::cellZoneSet::updateMesh(const mapPolyMesh& morphMap)
{
    // cellZone
    labelList newAddressing(addressing_.size());

    label n = 0;
    for (const label celli : addressing_)
    {
        label newCelli = morphMap.reverseCellMap()[celli];
        if (newCelli >= 0)
        {
            newAddressing[n] = newCelli;
            ++n;
        }
    }
    newAddressing.resize(n);

    addressing_.transfer(newAddressing);

    updateSet();
}


void Foam::cellZoneSet::writeDebug
(
    Ostream& os,
    const primitiveMesh& mesh,
    const label maxLen
) const
{
    cellSet::writeDebug(os, mesh, maxLen);
}


// ************************************************************************* //
