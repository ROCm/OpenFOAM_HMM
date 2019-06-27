/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2016 OpenFOAM Foundation
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
#include "bitSet.H"
#include "processorPolyPatch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ensightParts::ensightParts(const polyMesh& mesh)
:
    StorageType()
{
    recalculate(mesh);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ensightParts::recalculate(const polyMesh& mesh)
{
    StorageType::clear();

    label nPart = 0;

    // Track which cells are in a zone or not
    bitSet selection(mesh.nCells());

    // Do cell zones
    for (const cellZone& zn : mesh.cellZones())
    {
        if (zn.size())
        {
            selection.set(zn);
            this->append(new ensightPartCells(nPart++, mesh, zn));
        }
    }

    if (selection.none())
    {
        // No zones at all? - do entire mesh
        this->append(new ensightPartCells(nPart++, mesh));
    }
    else
    {
        // Flip from zoned to unzoned
        selection.flip();

        if (selection.any())
        {
            this->append(new ensightPartCells(nPart++, mesh, selection));
        }
    }


    // Do boundaries, skipping empty and processor patches
    for (const polyPatch& patch : mesh.boundaryMesh())
    {
        if (isA<processorPolyPatch>(patch))
        {
            // No processor patches
            break;
        }
        if (patch.size())
        {
            this->append(new ensightPartFaces(nPart++, mesh, patch));
        }
    }
}


void Foam::ensightParts::write(ensightGeoFile& os) const
{
    // Some feedback
    Info<< "Write geometry part (" << flush;

    for (const ensightPart& part : *this)
    {
        Info<< ' ' << part.index() << flush;
        part.write(os);
    }
    Info<< " )" << endl;
}


void Foam::ensightParts::writeSummary(Ostream& os) const
{
    for (const ensightPart& part : *this)
    {
        part.writeSummary(os);
    }
}


void Foam::ensightParts::dumpInfo(Ostream& os) const
{
    for (const ensightPart& part : *this)
    {
        part.dumpInfo(os);
    }
}


// ************************************************************************* //
