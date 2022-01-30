/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "rotatedBoxToCell.H"
#include "polyMesh.H"
#include "cellModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(rotatedBoxToCell, 0);
    addToRunTimeSelectionTable(topoSetSource, rotatedBoxToCell, word);
    addToRunTimeSelectionTable(topoSetSource, rotatedBoxToCell, istream);
    addToRunTimeSelectionTable(topoSetCellSource, rotatedBoxToCell, word);
    addToRunTimeSelectionTable(topoSetCellSource, rotatedBoxToCell, istream);
    addNamedToRunTimeSelectionTable
    (
        topoSetCellSource,
        rotatedBoxToCell,
        word,
        rotatedBox
    );
    addNamedToRunTimeSelectionTable
    (
        topoSetCellSource,
        rotatedBoxToCell,
        istream,
        rotatedBox
    );
}


Foam::topoSetSource::addToUsageTable Foam::rotatedBoxToCell::usage_
(
    rotatedBoxToCell::typeName,
    "\n    Usage: rotatedBoxToCell (originx originy originz)"
    " (ix iy iz) (jx jy jz) (kx ky kz)\n\n"
    "    Select all cells with cellCentre within parallelopiped\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::rotatedBoxToCell::combine(topoSet& set, const bool add) const
{
    // Define a cell for the box
    pointField boxPoints(8);
    boxPoints[0] = origin_;
    boxPoints[1] = origin_ + i_;
    boxPoints[2] = origin_ + i_ + j_;
    boxPoints[3] = origin_ + j_;
    boxPoints[4] = origin_ + k_;
    boxPoints[5] = origin_ + k_ + i_;
    boxPoints[6] = origin_ + k_ + i_ + j_;
    boxPoints[7] = origin_ + k_ + j_;

    labelList boxVerts(identity(8));

    const cellModel& hex = cellModel::ref(cellModel::HEX);

    // Get outwards pointing faces.
    faceList boxFaces(cellShape(hex, boxVerts).faces());

    // Precalculate normals
    vectorField boxFaceNormals(boxFaces.size());
    forAll(boxFaces, i)
    {
        boxFaceNormals[i] = boxFaces[i].areaNormal(boxPoints);

        //Pout<< "Face:" << i << " position:" << boxFaces[i].centre(boxPoints)
        //    << " normal:" << boxFaceNormals[i] << endl;
    }

    // Check whether cell centre is inside all faces of box.

    const pointField& ctrs = mesh_.cellCentres();

    forAll(ctrs, celli)
    {
        bool inside = true;

        forAll(boxFaces, i)
        {
            const face& f = boxFaces[i];

            if (((ctrs[celli] - boxPoints[f[0]]) & boxFaceNormals[i]) > 0)
            {
                inside = false;
                break;
            }
        }

        if (inside)
        {
            addOrDelete(set, celli, add);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rotatedBoxToCell::rotatedBoxToCell
(
    const polyMesh& mesh,
    const vector& origin,
    const vector& i,
    const vector& j,
    const vector& k
)
:
    topoSetCellSource(mesh),
    origin_(origin),
    i_(i),
    j_(j),
    k_(k)
{}


Foam::rotatedBoxToCell::rotatedBoxToCell
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    rotatedBoxToCell
    (
        mesh,
        dict.get<point>("origin"),
        dict.get<vector>("i"),
        dict.get<vector>("j"),
        dict.get<vector>("k")
    )
{}


Foam::rotatedBoxToCell::rotatedBoxToCell(const polyMesh& mesh, Istream& is)
:
    topoSetCellSource(mesh),
    origin_(is),
    i_(is),
    j_(is),
    k_(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::rotatedBoxToCell::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if (action == topoSetSource::ADD || action == topoSetSource::NEW)
    {
        if (verbose_)
        {
            Info<< "    Adding cells with centre within rotated box"
                << endl;
        }

        combine(set, true);
    }
    else if (action == topoSetSource::SUBTRACT)
    {
        if (verbose_)
        {
            Info<< "    Removing cells with centre within rotated box"
                << endl;
        }

        combine(set, false);
    }
}


// ************************************************************************* //
