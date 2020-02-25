/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "MeshedSurface.H"
#include "ListOps.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Face>
void Foam::MeshedSurface<Face>::checkZones()
{
    // extra safety, ensure we have at some zones
    // and they cover all the faces - fix start silently

    auto& zones = this->storedZones();

    label count = 0;
    for (surfZone& zn : zones)
    {
        zn.start() = count;
        count += zn.size();
    }

    if (!zones.empty())
    {
        if (count < this->size())
        {
            WarningInFunction
                << "more faces " << this->size() << " than zones " << count
                << " ... extending final zone" << nl;

            zones.last().size() += count - this->size();
        }
        else if (count > this->size())
        {
            FatalErrorInFunction
                << "more zones " << count << " than faces " << this->size()
                << exit(FatalError);
        }
    }
}


template<class Face>
void Foam::MeshedSurface<Face>::sortFacesAndStore
(
    DynamicList<Face>& unsortedFaces,
    DynamicList<label>& zoneIds,
    DynamicList<label>& elemIds,
    bool sorted
)
{
    // Basic sanity check
    const label nInputFaces = unsortedFaces.size();

    if (sorted || zoneIds.size() != nInputFaces)
    {
        // Sorting not required or not possible
        zoneIds.clear();
        sorted = true;
    }

    if (elemIds.size() != nInputFaces)
    {
        elemIds.clear();
    }

    if (sorted)
    {
        // No additional sorting required
        this->storedFaces().transfer(unsortedFaces);
        this->storedFaceIds().transfer(elemIds);
        return;
    }

    // The sorted order, based on zone-ids

    labelList faceMap;
    Foam::sortedOrder(zoneIds, faceMap);
    zoneIds.clear();

    auto& newFaces = this->storedFaces();
    newFaces.resize(nInputFaces);

    // Faces in sorted order
    forAll(newFaces, facei)
    {
        // Can use transfer, faceMap is unique
        newFaces[facei].transfer(unsortedFaces[faceMap[facei]]);
    }

    auto& newFaceIds = this->storedFaceIds();
    newFaceIds.resize(elemIds.size());

    // Element ids in sorted order
    forAll(newFaceIds, facei)
    {
        newFaceIds[facei] = elemIds[faceMap[facei]];
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
void Foam::MeshedSurface<Face>::addZones
(
    const UList<surfZone>& srfZones,
    const bool cullEmpty
)
{
    auto& zones = this->storedZones();
    zones.resize(zones.size());

    label nZone = 0;

    forAll(zones, zonei)
    {
        if (srfZones[zonei].size() || !cullEmpty)
        {
            zones[nZone] = surfZone(srfZones[zonei], nZone);
            ++nZone;
        }
    }

    zones.resize(nZone);
}


template<class Face>
void Foam::MeshedSurface<Face>::addZones
(
    const labelUList& sizes,
    const UList<word>& names,
    const bool cullEmpty
)
{
    auto& zones = this->storedZones();
    zones.resize(sizes.size());

    label start = 0;
    label nZone = 0;

    forAll(zones, zonei)
    {
        if (sizes[zonei] || !cullEmpty)
        {
            zones[nZone] = surfZone
            (
                names[zonei],
                sizes[zonei],
                start,
                nZone
            );
            start += sizes[zonei];
            ++nZone;
        }
    }

    zones.resize(nZone);
}


template<class Face>
void Foam::MeshedSurface<Face>::addZones
(
    const labelUList& sizes,
    const bool cullEmpty
)
{
    auto& zones = this->storedZones();
    zones.resize(sizes.size());

    label start = 0;
    label nZone = 0;

    forAll(zones, zonei)
    {
        if (sizes[zonei] || !cullEmpty)
        {
            zones[nZone] = surfZone
            (
                surfZone::defaultName(nZone),
                sizes[zonei],
                start,
                nZone
            );
            start += sizes[zonei];
            ++nZone;
        }
    }

    zones.resize(nZone);
}


template<class Face>
bool Foam::MeshedSurface<Face>::addZonesToFaces()
{
    // Normally a no-op, only the specializations are used.
    return false;
}


template<class Face>
void Foam::MeshedSurface<Face>::removeZones()
{
    this->storedZones().clear();
}


// ************************************************************************* //
