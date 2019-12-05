/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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

#include "primitiveMesh.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelListList& Foam::primitiveMesh::cellPoints() const
{
    if (!cpPtr_)
    {
        if (debug)
        {
            Pout<< "primitiveMesh::cellPoints() : "
                << "calculating cellPoints" << endl;

            if (debug == -1)
            {
                // For checking calls:abort so we can quickly hunt down
                // origin of call
                FatalErrorInFunction
                    << abort(FatalError);
            }
        }

        // Invert pointCells
        cpPtr_ = new labelListList(nCells());
        invertManyToMany(nCells(), pointCells(), *cpPtr_);
    }

    return *cpPtr_;
}


const Foam::labelList& Foam::primitiveMesh::cellPoints
(
    const label celli,
    labelHashSet& set,
    DynamicList<label>& storage
) const
{
    if (hasCellPoints())
    {
        return cellPoints()[celli];
    }

    const faceList& fcs = faces();
    const labelList& cFaces = cells()[celli];

    set.clear();

    for (const label facei : cFaces)
    {
        set.insert(fcs[facei]);
    }

    storage.clear();
    if (set.size() > storage.capacity())
    {
        storage.setCapacity(set.size());
    }

    for (const label pointi : set)
    {
        storage.append(pointi);
    }

    return storage;
}


const Foam::labelList& Foam::primitiveMesh::cellPoints(const label celli) const
{
    return cellPoints(celli, labelSet_, labels_);
}


// ************************************************************************* //
