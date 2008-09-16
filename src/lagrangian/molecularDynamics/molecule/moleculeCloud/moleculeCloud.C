/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*----------------------------------------------------------------------------*/

#include "moleculeCloud.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineParticleTypeNameAndDebug(molecule, 0);
    defineTemplateTypeNameAndDebug(Cloud<molecule>, 0);
};

Foam::scalar Foam::moleculeCloud::kb = 1.380650277e-23;

Foam::scalar Foam::moleculeCloud::elementaryCharge = 1.602176487e-19;

Foam::scalar Foam::moleculeCloud::vacuumPermittivity = 8.854187817e-12;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::moleculeCloud::moleculeCloud(const polyMesh& mesh)
:
    Cloud<molecule>(mesh, "moleculeCloud", false),
    mesh_(mesh),
    referredInteractionList_(*this)
    cellOccupancy_(mesh_.nCells())
{
    molecule::readFields(*this);

    buildCellOccupancy();

    removeHighEnergyOverlaps();

    calculateForce();
}


// Foam::moleculeCloud::moleculeCloud
// (
//     const polyMesh& mesh,
//     label nMol,
//     const labelField& id,
//     const scalarField& mass,
//     const vectorField& positions,
//     const labelField& cells,
//     const vectorField& U,
//     const vectorField& A,
//     const labelField& tethered,
//     const vectorField& tetherPositions
// )
// :
//     Cloud<molecule>(mesh, "moleculeCloud", false),
//     mesh_(mesh),
//     referredInteractionList_(*this)
// {
//     // Do I need to read the fields if I'm just about to clear them?
//     molecule::readFields(*this);

//     clear();

//     // This clear() is here for the moment to stop existing files
//     // being appended to, this would be better accomplished by getting
//     // mesh.removeFiles(mesh.instance()); (or equivalent) to work.

//     int i;

//     const Cloud<molecule>& cloud = *this;

//     for (i=0; i<nMol; i++)
//     {
//         addParticle
//         (
//             new molecule
//             (
//                 cloud,
//                 positions[i],
//                 cells[i],
//                 mass[i],
//                 U[i],
//                 A[i],
//                 tetherPositions[i],
//                 tethered[i],
//                 id[i]
//             )
//         );
//     }
// };


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::moleculeCloud::buildCellOccupancy()
{
    forAll(cellOccupancy_, cO)
    {
        cellOccupancy_[cO].clear();
    }

    iterator mol(this->begin());

    for
    (
        mol = this->begin();
        mol != this->end();
        ++mol
    )
    {
        cellOccupancy_[mol().cell()].append(&mol());
    }

    forAll(cellOccupancy_, cO)
    {
        cellOccupancy_[cO].shrink();
    }

    il_ril().referMolecules(cellOccupancy_);
}

void Foam::moleculeCloud::evolve()
{
    molecule::trackData td0(*this, 0);
    molecule::trackData td1(*this, 1);
    molecule::trackData td2(*this, 2);

    Cloud<molecule>::move(td0);

    calculateForce();

    Cloud<molecule>::move(td1);

    Cloud<molecule>::move(td2);

    Cloud<molecule>::move(td0);
}


void Foam::moleculeCloud::writeFields() const
{
    molecule::writeFields(*this);
}


// ************************************************************************* //
