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


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::moleculeCloud::buildConstProps()
{
    Info<< nl << "Reading moleculeProperties dictionary." << endl;

    const List<word>& idList(pot_.idList());

    constPropList_.setSize(idList.size());

    const List<word>& siteIdList(pot_.siteIdList());

    IOdictionary moleculePropertiesDict
    (
        IOobject
        (
            "moleculeProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    forAll(idList, i)
    {
        const word& id(idList[i]);

        const dictionary& molDict(moleculePropertiesDict.subDict(id));

        List<word> siteIdNames = molDict.lookup("siteIds");

        List<label> siteIds(siteIdNames.size());

        forAll(siteIdNames, sI)
        {
            const word& siteId = siteIdNames[sI];

            siteIds[sI] = findIndex(siteIdList, siteId);

            if (siteIds[sI] == -1)
            {
                FatalErrorIn("moleculeCloud.C") << nl
                    << siteId << " site not found."
                    << nl << abort(FatalError);
            }
        }

        molecule::constantProperties& constProp = constPropList_[i];

        constProp = molecule::constantProperties(molDict);

        constProp.siteIds() = siteIds;

        Info<< "sites " << constProp.siteIds() << endl;
    }
}


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

    il_.ril().referMolecules(cellOccupancy());
}


void Foam::moleculeCloud::calculatePairForce()
{
    iterator mol(this->begin());

    {
        // Real-Real interactions

        molecule* molI = &mol();

        molecule* molJ = &mol();

        const directInteractionList& dil(il_.dil());

        forAll(dil, d)
        {
            forAll(cellOccupancy_[d],cellIMols)
            {
                molI = cellOccupancy_[d][cellIMols];

                forAll(dil[d], interactingCells)
                {
                    List< molecule* > cellJ =
                    cellOccupancy_[dil[d][interactingCells]];

                    forAll(cellJ, cellJMols)
                    {
                        molJ = cellJ[cellJMols];

                        evaluatePair(molI, molJ);
                    }
                }

                forAll(cellOccupancy_[d],cellIOtherMols)
                {
                    molJ = cellOccupancy_[d][cellIOtherMols];

                    if (molJ > molI)
                    {
                        evaluatePair(molI, molJ);
                    }
                }
            }
        }
    }

    {
        // Real-Referred interactions

        molecule* molK = &mol();

        referredCellList& ril(il_.ril());

        forAll(ril, r)
        {
            const List<label>& realCells = ril[r].realCellsForInteraction();

            forAll(ril[r], refMols)
            {
                referredMolecule* molL(&(ril[r][refMols]));

                forAll(realCells, rC)
                {
                    List<molecule*> cellK = cellOccupancy_[realCells[rC]];

                    forAll(cellK, cellKMols)
                    {
                        molK = cellK[cellKMols];

                        evaluatePair(molK, molL);
                    }
                }
            }
        }
    }
}


void Foam::moleculeCloud::calculateElectrostaticForce()
{
    Info<< "Electrostatic forces" << endl;


}


void Foam::moleculeCloud::calculateTetherForce()
{
    const tetherPotentialList& tetherPot(pot_.tetherPotentials());

    iterator mol(this->begin());

    for (mol = this->begin(); mol != this->end(); ++mol)
    {
        if (mol().tethered())
        {
            vector rIT = mol().position() - mol().specialPosition();

            label idI = mol().id();

            scalar massI = constProps(idI).mass();

            vector fIT = tetherPot.force(idI, rIT);

            mol().a() += fIT/massI;

            mol().potentialEnergy() += tetherPot.energy(idI, rIT);

            mol().rf() += rIT*fIT;
        }
    }
}


void Foam::moleculeCloud::calculateExternalForce()
{
    iterator mol(this->begin());

    // Info<< "Warning! Includes dissipation term!" << endl;

    for (mol = this->begin(); mol != this->end(); ++mol)
    {
        mol().a() += pot_.gravity();

        // mol().a() += -1.0 * mol().U() /constProps(mol.id()).mass()
    }
}


void Foam::moleculeCloud::removeHighEnergyOverlaps()
{
    Info << nl << "Removing high energy overlaps, removal order:";

    forAll(pot_.removalOrder(), rO)
    {
        Info << " " << pot_.idList()[pot_.removalOrder()[rO]];
    }

    Info << nl ;

    label initialSize = this->size();

    buildCellOccupancy();

    // Real-Real interaction
    iterator mol(this->begin());

    {
        mol = this->begin();

        molecule* molI = &mol();

        molecule* molJ = &mol();

        DynamicList<molecule*> molsToDelete;

        const directInteractionList& dil(il_.dil());

        forAll(dil, d)
        {
            forAll(cellOccupancy_[d],cellIMols)
            {
                molI = cellOccupancy_[d][cellIMols];

                forAll(dil[d], interactingCells)
                {
                    List< molecule* > cellJ =
                    cellOccupancy_[dil[d][interactingCells]];

                    forAll(cellJ, cellJMols)
                    {
                        molJ = cellJ[cellJMols];

                        if (evaluatePotentialLimit(molI, molJ))
                        {
                            label idI = molI->id();

                            label idJ = molJ->id();

                            if
                            (
                                idI == idJ
                             || findIndex(pot_.removalOrder(), idJ)
                                    < findIndex(pot_.removalOrder(), idI)
                            )
                            {
                                if (findIndex(molsToDelete, molJ) == -1)
                                {
                                    molsToDelete.append(molJ);
                                }
                            }
                            else if (findIndex(molsToDelete, molI) == -1)
                            {
                                molsToDelete.append(molI);
                            }
                        }
                    }
                }
            }

            forAll(cellOccupancy_[d],cellIOtherMols)
            {
                molJ = cellOccupancy_[d][cellIOtherMols];

                if (molJ > molI)
                {
                    if (evaluatePotentialLimit(molI, molJ))
                    {
                        label idI = molI->id();

                        label idJ = molJ->id();

                        if
                        (
                            idI == idJ
                         || findIndex(pot_.removalOrder(), idJ)
                                < findIndex(pot_.removalOrder(), idI)
                        )
                        {
                            if (findIndex(molsToDelete, molJ) == -1)
                            {
                                molsToDelete.append(molJ);
                            }
                        }
                        else if (findIndex(molsToDelete, molI) == -1)
                        {
                            molsToDelete.append(molI);
                        }
                    }
                }
            }
        }

        forAll (molsToDelete, mTD)
        {
            deleteParticle(*(molsToDelete[mTD]));
        }
    }

    buildCellOccupancy();

    // Real-Referred interaction

    {
        molecule* molK = &mol();

        DynamicList<molecule*> molsToDelete;

        referredCellList& ril(il_.ril());

        forAll(ril, r)
        {
            referredCell& refCell = ril[r];

            forAll(refCell, refMols)
            {
                referredMolecule* molL = &(refCell[refMols]);

                List <label> realCells = refCell.realCellsForInteraction();

                forAll(realCells, rC)
                {
                    label cellK = realCells[rC];

                    List<molecule*> cellKMols = cellOccupancy_[cellK];

                    forAll(cellKMols, cKM)
                    {
                        molK = cellKMols[cKM];

                        if (evaluatePotentialLimit(molK, molL))
                        {
                            label idK = molK->id();

                            label idL = molL->id();

                            if
                            (
                                findIndex(pot_.removalOrder(), idK)
                                    < findIndex(pot_.removalOrder(), idL)
                            )
                            {
                                if (findIndex(molsToDelete, molK) == -1)
                                {
                                    molsToDelete.append(molK);
                                }
                            }

                            else if
                            (
                                findIndex(pot_.removalOrder(), idK)
                                    == findIndex(pot_.removalOrder(), idL)
                            )
                            {
                                if
                                (
                                    Pstream::myProcNo() == refCell.sourceProc()
                                 && cellK <= refCell.sourceCell()
                                )
                                {
                                    if (findIndex(molsToDelete, molK) == -1)
                                    {
                                        molsToDelete.append(molK);
                                    }
                                }

                                else if
                                (
                                    Pstream::myProcNo() < refCell.sourceProc()
                                )
                                {
                                    if (findIndex(molsToDelete, molK) == -1)
                                    {
                                        molsToDelete.append(molK);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        forAll (molsToDelete, mTD)
        {
            deleteParticle(*(molsToDelete[mTD]));
        }
    }

    buildCellOccupancy();

    label molsRemoved = initialSize - this->size();

    if (Pstream::parRun())
    {
        reduce(molsRemoved, sumOp<label>());
    }

    Info << tab << molsRemoved << " molecules removed" << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::moleculeCloud::moleculeCloud
(
    const polyMesh& mesh,
    const potential& pot
)
:
    Cloud<molecule>(mesh, "moleculeCloud", false),
    mesh_(mesh),
    pot_(pot),
    cellOccupancy_(mesh_.nCells()),
    il_(mesh_, pot_.pairPotentials().rCutMaxSqr(), false),
    constPropList_()
{
    molecule::readFields(*this);

    buildConstProps();

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

void Foam::moleculeCloud::evolve()
{
    molecule::trackData td0(*this, 0);
    Cloud<molecule>::move(td0);

    molecule::trackData td1(*this, 1);
    Cloud<molecule>::move(td1);

    molecule::trackData td2(*this, 2);
    Cloud<molecule>::move(td2);

    calculateForce();

    molecule::trackData td3(*this, 3);
    Cloud<molecule>::move(td3);
}


void Foam::moleculeCloud::calculateForce()
{
    buildCellOccupancy();

    iterator mol(this->begin());

    // Set accumulated quantities to zero
    for (mol = this->begin(); mol != this->end(); ++mol)
    {
        mol().siteForces() = vector::zero;

        mol().potentialEnergy() = 0.0;

        mol().rf() = tensor::zero;
    }

    calculatePairForce();

    calculateElectrostaticForce();

    calculateTetherForce();

    calculateExternalForce();
}


void Foam::moleculeCloud::applyConstraintsAndThermostats
(
    const scalar targetTemperature,
    const scalar measuredTemperature
)
{
    scalar temperatureCorrectionFactor =
        sqrt(targetTemperature/measuredTemperature);

    Info<< "----------------------------------------" << nl
        << "Temperature equilibration" << nl
        << "Target temperature = "
        << targetTemperature << nl
        << "Measured temperature = "
        << measuredTemperature << nl
        << "Temperature correction factor = "
        << temperatureCorrectionFactor << nl
        << "----------------------------------------"
        << endl;

    iterator mol(this->begin());

    for (mol = this->begin(); mol != this->end(); ++mol)
    {
        mol().v() *= temperatureCorrectionFactor;
    }
}


void Foam::moleculeCloud::writeFields() const
{
    molecule::writeFields(*this);
}


// ************************************************************************* //
