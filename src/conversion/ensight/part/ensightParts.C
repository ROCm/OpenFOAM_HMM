/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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

#include "ensightParts.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ensightParts::ensightParts(const polyMesh& mesh)
:
    StorageType()
{
    recalculate(mesh);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ensightParts::~ensightParts()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ensightParts::recalculate(const polyMesh& mesh)
{
    StorageType::clear();

    label nPart = 0;
    label nZoneCells = 0;

    // do cell zones
    forAll(mesh.cellZones(), zoneI)
    {
        const cellZone& cZone = mesh.cellZones()[zoneI];
        nZoneCells += cZone.size();

        if (cZone.size())
        {
            this->append(new ensightPartCells(nPart++, mesh, cZone));
        }
    }

    // collect unzoned cells

    // special case: no zones at all - do entire mesh
    if (nZoneCells == 0)
    {
        this->append(new ensightPartCells(nPart++, mesh));
    }
    else if (mesh.nCells() > nZoneCells)
    {
        // determine which cells are not in a cellZone
        labelList unzoned(mesh.nCells(), -1);

        forAll(mesh.cellZones(), zoneI)
        {
            const labelUList& idList = mesh.cellZones()[zoneI];

            forAll(idList, i)
            {
                unzoned[idList[i]] = idList[i];
            }
        }

        label nUnzoned = 0;
        forAll(unzoned, i)
        {
            if (unzoned[i] < 0)
            {
                unzoned[nUnzoned] = i;
                nUnzoned++;
            }
        }
        unzoned.setSize(nUnzoned);

        if (unzoned.size())
        {
            this->append(new ensightPartCells(nPart++, mesh, unzoned));
        }
    }


    // do boundaries, skipping empty and processor patches
    forAll(mesh.boundaryMesh(), patchi)
    {
        const polyPatch& patch = mesh.boundaryMesh()[patchi];
        if (patch.size() && !isA<processorPolyPatch>(patch))
        {
            this->append(new ensightPartFaces(nPart++, mesh, patch));
        }
    }
}


void Foam::ensightParts::write(ensightGeoFile& os) const
{
    // Some feedback
    Info<< "Write geometry part (" << flush;

    forAllConstIter(StorageType, *this, iter)
    {
        Info<< ' ' << (*iter).index() << flush;
        (*iter).write(os);
    }
    Info<< " )" << endl;
}


void Foam::ensightParts::writeSummary(Ostream& os) const
{
    forAllConstIter(StorageType, *this, iter)
    {
        (*iter).writeSummary(os);
    }
}


void Foam::ensightParts::dumpInfo(Ostream& os) const
{
    forAllConstIter(StorageType, *this, iter)
    {
        (*iter).dumpInfo(os);
    }
}


// ************************************************************************* //
