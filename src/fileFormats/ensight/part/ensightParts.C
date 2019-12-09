/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2019 OpenCFD Ltd.
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
#include "emptyPolyPatch.H"
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

    // Do all cell zones
    for (const cellZone& zn : mesh.cellZones())
    {
        if (returnReduce(!zn.empty(), orOp<bool>()))
        {
            selection.set(zn);
            this->append(new ensightPartCells(nPart++, mesh, zn));
        }
    }

    if (!nPart)
    {
        // No zones at all? - do entire mesh. Name as per ensightMesh
        this->append(new ensightPartCells(nPart++, mesh, "internalMesh"));
    }
    else
    {
        // Flip from zoned to unzoned
        selection.flip();

        if (returnReduce(selection.any(), orOp<bool>()))
        {
            this->append
            (
                new ensightPartCells(nPart++, mesh, selection, "__internal__")
            );
        }
    }


    // Do boundaries, skipping empty and processor patches
    for (const polyPatch& p : mesh.boundaryMesh())
    {
        if (isA<processorPolyPatch>(p))
        {
            // No processor patches
            break;
        }

        if (returnReduce(!p.empty(), orOp<bool>()))
        {
            this->append(new ensightPartFaces(nPart++, mesh, p));
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
