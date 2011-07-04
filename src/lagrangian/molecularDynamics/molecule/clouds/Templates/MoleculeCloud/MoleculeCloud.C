/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2011 OpenCFD Ltd.
     \\/     M anipulation  |
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

\*----------------------------------------------------------------------------*/

#include "MoleculeCloud.H"
#include "fvMesh.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class MoleculeType>
void Foam::MoleculeCloud<MoleculeType>::buildConstProps()
{
    const List<word>& idList(pot_.idList());

    constPropList_.setSize(idList.size());

    const List<word>& siteIdList(pot_.siteIdList());

    IOdictionary molPropertiesDict
    (
        IOobject
        (
            this->name() + "Properties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    );

    Info<< nl << "Reading " << molPropertiesDict.name() << endl;

    forAll(idList, i)
    {
        const word& id = idList[i];

        const dictionary& molDict(molPropertiesDict.subDict(id));

        List<word> siteIdNames = molDict.lookup("siteIds");

        List<label> siteIds(siteIdNames.size());

        forAll(siteIdNames, sI)
        {
            const word& siteId = siteIdNames[sI];

            siteIds[sI] = findIndex(siteIdList, siteId);

            if (siteIds[sI] == -1)
            {
                FatalErrorIn
                (
                    "void Foam::MoleculeCloud<MoleculeType>::buildConstProps()"
                )
                    << siteId << " site not found."
                    << nl << abort(FatalError);
            }
        }

        constPropList_[i] = typename MoleculeType::constantProperties
        (
            molDict,
            siteIds
        );
    }
}


template<class MoleculeType>
void Foam::MoleculeCloud<MoleculeType>::setSiteSizesAndPositions()
{
    forAllIter(typename MoleculeCloud<MoleculeType>, *this, mol)
    {
        const typename MoleculeType::constantProperties& cP
        (
            constProps(mol().id())
        );

        mol().setSiteSizes(cP.nSites());

        mol().setSitePositions(cP);
    }
}


template<class MoleculeType>
void Foam::MoleculeCloud<MoleculeType>::buildCellOccupancy()
{
    forAll(cellOccupancy_, cO)
    {
        cellOccupancy_[cO].clear();
    }

    forAllIter(typename MoleculeCloud<MoleculeType>, *this, mol)
    {
        cellOccupancy_[mol().cell()].append(&mol());
    }

    forAll(cellOccupancy_, cO)
    {
        cellOccupancy_[cO].shrink();
    }
}


template<class MoleculeType>
void Foam::MoleculeCloud<MoleculeType>::calculatePairForce()
{
    PstreamBuffers pBufs(Pstream::nonBlocking);

    // Start sending referred data
    il_.sendReferredData(cellOccupancy(), pBufs);

    MoleculeType* molI = NULL;
    MoleculeType* molJ = NULL;

    {
        // Real-Real interactions

        const labelListList& dil = il_.dil();

        forAll(dil, d)
        {
            forAll(cellOccupancy_[d],cellIMols)
            {
                molI = cellOccupancy_[d][cellIMols];

                forAll(dil[d], interactingCells)
                {
                    List<MoleculeType*> cellJ =
                        cellOccupancy_[dil[d][interactingCells]];

                    forAll(cellJ, cellJMols)
                    {
                        molJ = cellJ[cellJMols];

                        evaluatePair(*molI, *molJ);
                    }
                }

                forAll(cellOccupancy_[d], cellIOtherMols)
                {
                    molJ = cellOccupancy_[d][cellIOtherMols];

                    if (molJ > molI)
                    {
                        evaluatePair(*molI, *molJ);
                    }
                }
            }
        }
    }

    // Receive referred data
    il_.receiveReferredData(pBufs);

    {
        // Real-Referred interactions

        const labelListList& ril = il_.ril();

        List<IDLList<MoleculeType> >& referredMols = il_.referredParticles();

        forAll(ril, r)
        {
            const List<label>& realCells = ril[r];

            IDLList<MoleculeType>& refMols = referredMols[r];

            forAllIter
            (
                typename IDLList<MoleculeType>,
                refMols,
                refMol
            )
            {
                forAll(realCells, rC)
                {
                    List<MoleculeType*> cellI = cellOccupancy_[realCells[rC]];

                    forAll(cellI, cellIMols)
                    {
                        molI = cellI[cellIMols];

                        evaluatePair(*molI, refMol());
                    }
                }
            }
        }
    }
}


template<class MoleculeType>
void Foam::MoleculeCloud<MoleculeType>::calculateTetherForce()
{
    const tetherPotentialList& tetherPot(pot_.tetherPotentials());

    forAllIter(typename MoleculeCloud<MoleculeType>, *this, mol)
    {
        if (mol().tethered())
        {
            vector rIT = mol().position() - mol().specialPosition();

            label idI = mol().id();

            scalar massI = constProps(idI).mass();

            vector fIT = tetherPot.force(idI, rIT);

            mol().a() += fIT/massI;

            mol().potentialEnergy() += tetherPot.energy(idI, rIT);

            // What to do here?
            // mol().rf() += rIT*fIT;
        }
    }
}


template<class MoleculeType>
void Foam::MoleculeCloud<MoleculeType>::calculateExternalForce()
{
    forAllIter(typename MoleculeCloud<MoleculeType>, *this, mol)
    {
        mol().a() += pot_.gravity();
    }
}


template<class MoleculeType>
void Foam::MoleculeCloud<MoleculeType>::removeHighEnergyOverlaps()
{
    Info<< nl << "Removing high energy overlaps, limit = "
        << pot_.potentialEnergyLimit()
        << nl << "Removal order:";

    forAll(pot_.removalOrder(), rO)
    {
        Info<< ' ' << pot_.idList()[pot_.removalOrder()[rO]];
    }

    Info<< nl ;

    label initialSize = this->size();

    buildCellOccupancy();

    // Real-Real interaction

    MoleculeType* molI = NULL;
    MoleculeType* molJ = NULL;

    {
        DynamicList<MoleculeType*> molsToDelete;

        const labelListList& dil(il_.dil());

        forAll(dil, d)
        {
            forAll(cellOccupancy_[d], cellIMols)
            {
                molI = cellOccupancy_[d][cellIMols];

                forAll(dil[d], interactingCells)
                {
                    List<MoleculeType*> cellJ =
                        cellOccupancy_[dil[d][interactingCells]];

                    forAll(cellJ, cellJMols)
                    {
                        molJ = cellJ[cellJMols];

                        if (evaluatePotentialLimit(*molI, *molJ))
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

                forAll(cellOccupancy_[d], cellIOtherMols)
                {
                    molJ = cellOccupancy_[d][cellIOtherMols];

                    if (molJ > molI)
                    {
                        if (evaluatePotentialLimit(*molI, *molJ))
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
        }

        forAll(molsToDelete, mTD)
        {
            deleteParticle(*(molsToDelete[mTD]));
        }
    }

    buildCellOccupancy();

    PstreamBuffers pBufs(Pstream::nonBlocking);

    // Start sending referred data
    il_.sendReferredData(cellOccupancy(), pBufs);

        // Receive referred data
    il_.receiveReferredData(pBufs);

    // Real-Referred interaction

    {
        DynamicList<MoleculeType*> molsToDelete;

        const labelListList& ril(il_.ril());

        List<IDLList<MoleculeType> >& referredMols = il_.referredParticles();

        forAll(ril, r)
        {
            IDLList<MoleculeType>& refMols = referredMols[r];

            forAllIter
            (
                typename IDLList<MoleculeType>,
                refMols,
                refMol
            )
            {
                molJ = &refMol();

                const List<label>& realCells = ril[r];

                forAll(realCells, rC)
                {
                    label cellI = realCells[rC];

                    List<MoleculeType*> cellIMols = cellOccupancy_[cellI];

                    forAll(cellIMols, cIM)
                    {
                        molI = cellIMols[cIM];

                        if (evaluatePotentialLimit(*molI, *molJ))
                        {
                            label idI = molI->id();

                            label idJ = molJ->id();

                            if
                            (
                                findIndex(pot_.removalOrder(), idI)
                              < findIndex(pot_.removalOrder(), idJ)
                            )
                            {
                                if (findIndex(molsToDelete, molI) == -1)
                                {
                                    molsToDelete.append(molI);
                                }
                            }
                            else if
                            (
                                findIndex(pot_.removalOrder(), idI)
                             == findIndex(pot_.removalOrder(), idJ)
                            )
                            {
                                // Remove one of the molecules
                                // arbitrarily, assuring that a
                                // consistent decision is made for
                                // both real-referred pairs.

                                if (molI->origId() > molJ->origId())
                                {
                                    if (findIndex(molsToDelete, molI) == -1)
                                    {
                                        molsToDelete.append(molI);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        forAll(molsToDelete, mTD)
        {
            deleteParticle(*(molsToDelete[mTD]));
        }
    }

    buildCellOccupancy();

    // Start sending referred data
    il_.sendReferredData(cellOccupancy(), pBufs);

    // Receive referred data
    il_.receiveReferredData(pBufs);

    label molsRemoved = initialSize - this->size();

    if (Pstream::parRun())
    {
        reduce(molsRemoved, sumOp<label>());
    }

    Info<< tab << molsRemoved << " molecules removed" << endl;
}


template<class MoleculeType>
void Foam::MoleculeCloud<MoleculeType>::initialiseMolecules
(
    const IOdictionary& mdInitialiseDict
)
{
    Info<< nl
        << "Initialising molecules in each zone specified in "
        << mdInitialiseDict.name()
        << endl;

    const cellZoneMesh& cellZones = mesh_.cellZones();

    if (!cellZones.size())
    {
        FatalErrorIn
        (
            "void Foam::MoleculeCloud<MoleculeType>::initialiseMolecules"
        )
            << "No cellZones found in the mesh."
            << abort(FatalError);
    }

    forAll(cellZones, z)
    {
        const cellZone& zone(cellZones[z]);

        if (zone.size())
        {
            if (!mdInitialiseDict.found(zone.name()))
            {
                Info<< "No specification subDictionary for zone "
                    << zone.name() << " found, skipping." << endl;
            }
            else
            {
                const dictionary& zoneDict =
                    mdInitialiseDict.subDict(zone.name());

                const scalar temperature
                (
                    readScalar(zoneDict.lookup("temperature"))
                );

                const vector bulkVelocity(zoneDict.lookup("bulkVelocity"));

                List<word> latticeIds
                (
                    zoneDict.lookup("latticeIds")
                );

                List<vector> latticePositions
                (
                    zoneDict.lookup("latticePositions")
                );

                if (latticeIds.size() != latticePositions.size())
                {
                    FatalErrorIn
                    (
                        "Foam::MoleculeCloud<MoleculeType>::initialiseMolecules"
                    )
                        << "latticeIds and latticePositions must be the same "
                        << " size." << nl
                        << abort(FatalError);
                }

                diagTensor latticeCellShape
                (
                    zoneDict.lookup("latticeCellShape")
                );

                scalar latticeCellScale = 0.0;

                if (zoneDict.found("numberDensity"))
                {
                    scalar numberDensity = readScalar
                    (
                        zoneDict.lookup("numberDensity")
                    );

                    if (numberDensity < VSMALL)
                    {
                        WarningIn
                        (
                            "MoleculeCloud<MoleculeType>::initialiseMolecules"
                        )
                            << "numberDensity too small, not filling zone "
                            << zone.name() << endl;

                        continue;
                    }

                    latticeCellScale = pow
                    (
                        latticeIds.size()/(det(latticeCellShape)*numberDensity),
                        (1.0/3.0)
                    );
                }
                else if (zoneDict.found("massDensity"))
                {
                    scalar unitCellMass = 0.0;

                    forAll(latticeIds, i)
                    {
                        label id = findIndex(pot_.idList(), latticeIds[i]);

                        const typename MoleculeType::constantProperties& cP
                        (
                            constProps(id)
                        );

                        unitCellMass += cP.mass();
                    }

                    scalar massDensity = readScalar
                    (
                        zoneDict.lookup("massDensity")
                    );

                    if (massDensity < VSMALL)
                    {
                        WarningIn
                        (
                            "MoleculeCloud<MoleculeType>::initialiseMolecules"
                        )
                            << "massDensity too small, not filling zone "
                            << zone.name() << endl;

                        continue;
                    }


                    latticeCellScale = pow
                    (
                        unitCellMass/(det(latticeCellShape)*massDensity),
                        (1.0/3.0)
                    );
                }
                else
                {
                    FatalErrorIn
                    (
                        "Foam::MoleculeCloud<MoleculeType>::initialiseMolecules"
                    )
                        << "massDensity or numberDensity not specified " << nl
                        << abort(FatalError);
                }

                latticeCellShape *= latticeCellScale;

                vector anchor(zoneDict.lookup("anchor"));

                bool tethered = false;

                if (zoneDict.found("tetherSiteIds"))
                {
                    tethered = bool
                    (
                        List<word>(zoneDict.lookup("tetherSiteIds")).size()
                    );
                }

                const vector orientationAngles
                (
                    zoneDict.lookup("orientationAngles")
                );

                scalar phi(orientationAngles.x()*pi/180.0);

                scalar theta(orientationAngles.y()*pi/180.0);

                scalar psi(orientationAngles.z()*pi/180.0);

                const tensor R
                (
                    cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi),
                    cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi),
                    sin(psi)*sin(theta),
                    - sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi),
                    - sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi),
                    cos(psi)*sin(theta),
                    sin(theta)*sin(phi),
                    - sin(theta)*cos(phi),
                    cos(theta)
                );

                // Find the optimal anchor position.  Finding the approximate
                // mid-point of the zone of cells and snapping to the nearest
                // lattice location.

                vector zoneMin = VGREAT*vector::one;

                vector zoneMax = -VGREAT*vector::one;

                forAll(zone, cell)
                {
                    const point cellCentre = mesh_.cellCentres()[zone[cell]];

                    if (cellCentre.x() > zoneMax.x())
                    {
                        zoneMax.x() = cellCentre.x();
                    }
                    if (cellCentre.x() < zoneMin.x())
                    {
                        zoneMin.x() = cellCentre.x();
                    }
                    if (cellCentre.y() > zoneMax.y())
                    {
                        zoneMax.y() = cellCentre.y();
                    }
                    if (cellCentre.y() < zoneMin.y())
                    {
                        zoneMin.y() = cellCentre.y();
                    }
                    if (cellCentre.z() > zoneMax.z())
                    {
                        zoneMax.z() = cellCentre.z();
                    }
                    if (cellCentre.z() < zoneMin.z())
                    {
                        zoneMin.z() = cellCentre.z();
                    }
                }

                point zoneMid = 0.5*(zoneMin + zoneMax);

                point latticeMid = inv(latticeCellShape) & (R.T()
                & (zoneMid - anchor));

                point latticeAnchor
                (
                    label(latticeMid.x() + 0.5*sign(latticeMid.x())),
                    label(latticeMid.y() + 0.5*sign(latticeMid.y())),
                    label(latticeMid.z() + 0.5*sign(latticeMid.z()))
                );

                anchor += (R & (latticeCellShape & latticeAnchor));

                // Continue trying to place molecule as long as at
                // least one molecule is placed in each iteration.
                // The "|| totalZoneMols == 0" condition means that the
                // algorithm will continue if the origin is outside the
                // zone.

                label n = 0;

                label totalZoneMols = 0;

                label molsPlacedThisIteration = 0;

                while
                (
                    molsPlacedThisIteration != 0
                 || totalZoneMols == 0
                )
                {
                    label sizeBeforeIteration = this->size();

                    bool partOfLayerInBounds = false;

                    if (n == 0)
                    {
                        // Special treatment is required for the first position,
                        // i.e. iteration zero.

                        labelVector unitCellLatticePosition(0,0,0);

                        forAll(latticePositions, p)
                        {
                            label id = findIndex(pot_.idList(), latticeIds[p]);

                            const vector& latticePosition =
                                vector
                                (
                                    unitCellLatticePosition.x(),
                                    unitCellLatticePosition.y(),
                                    unitCellLatticePosition.z()
                                )
                              + latticePositions[p];

                            point globalPosition =
                                anchor
                              + (R & (latticeCellShape & latticePosition));

                            partOfLayerInBounds = mesh_.bounds().contains
                            (
                                globalPosition
                            );

                            label cell = -1;
                            label tetFace = -1;
                            label tetPt = -1;

                            mesh_.findCellFacePt
                            (
                                globalPosition,
                                cell,
                                tetFace,
                                tetPt
                            );

                            if (findIndex(zone, cell) != -1)
                            {
                                createMolecule
                                (
                                    globalPosition,
                                    cell,
                                    tetFace,
                                    tetPt,
                                    id,
                                    tethered,
                                    temperature,
                                    bulkVelocity
                                );
                            }
                        }
                    }
                    else
                    {
                        // Place top and bottom caps.

                        labelVector unitCellLatticePosition(0,0,0);

                        for
                        (
                            unitCellLatticePosition.z() = -n;
                            unitCellLatticePosition.z() <= n;
                            unitCellLatticePosition.z() += 2*n
                        )
                        {
                            for
                            (
                                unitCellLatticePosition.y() = -n;
                                unitCellLatticePosition.y() <= n;
                                unitCellLatticePosition.y()++
                            )
                            {
                                for
                                (
                                    unitCellLatticePosition.x() = -n;
                                    unitCellLatticePosition.x() <= n;
                                    unitCellLatticePosition.x()++
                                )
                                {
                                    forAll(latticePositions, p)
                                    {
                                        label id = findIndex
                                        (
                                            pot_.idList(),
                                            latticeIds[p]
                                        );

                                        const vector& latticePosition =
                                            vector
                                            (
                                                unitCellLatticePosition.x(),
                                                unitCellLatticePosition.y(),
                                                unitCellLatticePosition.z()
                                            )
                                          + latticePositions[p];

                                        point globalPosition =
                                            anchor
                                          + (
                                                R
                                              & (
                                                    latticeCellShape
                                                  & latticePosition
                                                )
                                            );

                                        partOfLayerInBounds =
                                            mesh_.bounds().contains
                                            (
                                                globalPosition
                                            );

                                        label cell = -1;
                                        label tetFace = -1;
                                        label tetPt = -1;

                                        mesh_.findCellFacePt
                                        (
                                            globalPosition,
                                            cell,
                                            tetFace,
                                            tetPt
                                        );

                                        if (findIndex(zone, cell) != -1)
                                        {
                                            createMolecule
                                            (
                                                globalPosition,
                                                cell,
                                                tetFace,
                                                tetPt,
                                                id,
                                                tethered,
                                                temperature,
                                                bulkVelocity
                                            );
                                        }
                                    }
                                }
                            }
                        }

                        for
                        (
                            unitCellLatticePosition.z() = -(n-1);
                            unitCellLatticePosition.z() <= (n-1);
                            unitCellLatticePosition.z()++
                        )
                        {
                            for (label iR = 0; iR <= 2*n -1; iR++)
                            {
                                unitCellLatticePosition.x() = n;

                                unitCellLatticePosition.y() = -n + (iR + 1);

                                for (label iK = 0; iK < 4; iK++)
                                {
                                    forAll(latticePositions, p)
                                    {
                                        label id = findIndex
                                        (
                                            pot_.idList(),
                                            latticeIds[p]
                                        );

                                        const vector& latticePosition =
                                            vector
                                            (
                                                unitCellLatticePosition.x(),
                                                unitCellLatticePosition.y(),
                                                unitCellLatticePosition.z()
                                            )
                                          + latticePositions[p];

                                        point globalPosition =
                                            anchor
                                          + (
                                                R
                                              & (
                                                   latticeCellShape
                                                 & latticePosition
                                                )
                                            );

                                        partOfLayerInBounds =
                                            mesh_.bounds().contains
                                            (
                                                globalPosition
                                            );

                                        label cell = -1;
                                        label tetFace = -1;
                                        label tetPt = -1;

                                        mesh_.findCellFacePt
                                        (
                                            globalPosition,
                                            cell,
                                            tetFace,
                                            tetPt
                                        );

                                        if (findIndex(zone, cell) != -1)
                                        {
                                            createMolecule
                                            (
                                                globalPosition,
                                                cell,
                                                tetFace,
                                                tetPt,
                                                id,
                                                tethered,
                                                temperature,
                                                bulkVelocity
                                            );
                                        }
                                    }

                                    unitCellLatticePosition =
                                        labelVector
                                        (
                                          - unitCellLatticePosition.y(),
                                            unitCellLatticePosition.x(),
                                            unitCellLatticePosition.z()
                                        );
                                }
                            }
                        }
                    }

                    if
                    (
                        totalZoneMols == 0
                     && !partOfLayerInBounds
                    )
                    {
                        WarningIn
                        (
                            "Foam::MoleculeCloud<MoleculeType>::"
                            "initialiseMolecules()"
                        )
                            << "A whole layer of unit cells was placed "
                            << "outside the bounds of the mesh, but no "
                            << "molecules have been placed in zone '"
                            << zone.name()
                            << "'.  This is likely to be because the zone "
                            << "has few cells ("
                            << zone.size()
                            << " in this case) and no lattice position "
                            << "fell inside them.  "
                            << "Aborting filling this zone."
                            << endl;

                        break;
                    }

                    molsPlacedThisIteration =
                        this->size() - sizeBeforeIteration;

                    totalZoneMols += molsPlacedThisIteration;

                    n++;
                }
            }
        }
    }
}


template<class MoleculeType>
void Foam::MoleculeCloud<MoleculeType>::createMolecule
(
    const point& position,
    label cell,
    label tetFace,
    label tetPt,
    label id,
    bool tethered,
    scalar temperature,
    const vector& bulkVelocity
)
{
    if (cell == -1)
    {
        mesh_.findCellFacePt(position, cell, tetFace, tetPt);
    }

    if (cell == -1)
    {
        FatalErrorIn("Foam::MoleculeCloud<MoleculeType>::createMolecule")
            << "Position specified does not correspond to a mesh cell." << nl
            << abort(FatalError);
    }

    point specialPosition(vector::zero);

    label special = 0;

    if (tethered)
    {
        specialPosition = position;

        special = MoleculeType::SPECIAL_TETHERED;
    }

    const typename MoleculeType::constantProperties& cP(constProps(id));

    vector v = equipartitionLinearVelocity(temperature, cP.mass());

    v += bulkVelocity;

    vector pi = vector::zero;

    tensor Q = I;

    if (!cP.pointMolecule())
    {
        pi = equipartitionAngularMomentum(temperature, cP);

        scalar phi(rndGen_.scalar01()*twoPi);

        scalar theta(rndGen_.scalar01()*twoPi);

        scalar psi(rndGen_.scalar01()*twoPi);

        Q = tensor
        (
            cos(psi)*cos(phi) - cos(theta)*sin(phi)*sin(psi),
            cos(psi)*sin(phi) + cos(theta)*cos(phi)*sin(psi),
            sin(psi)*sin(theta),
            - sin(psi)*cos(phi) - cos(theta)*sin(phi)*cos(psi),
            - sin(psi)*sin(phi) + cos(theta)*cos(phi)*cos(psi),
            cos(psi)*sin(theta),
            sin(theta)*sin(phi),
            - sin(theta)*cos(phi),
            cos(theta)
        );
    }

    addParticle
    (
        new MoleculeType
        (
            mesh_,
            position,
            cell,
            tetFace,
            tetPt,
            Q,
            v,
            vector::zero,
            pi,
            vector::zero,
            specialPosition,
            constProps(id),
            special,
            id
        )
    );
}


template<class MoleculeType>
Foam::label Foam::MoleculeCloud<MoleculeType>::nSites() const
{
    label n = 0;

    forAllConstIter(typename MoleculeCloud<MoleculeType>, *this, mol)
    {
        n += constProps(mol().id()).nSites();
    }

    return n;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class MoleculeType>
Foam::MoleculeCloud<MoleculeType>::MoleculeCloud
(
    const word& cloudName,
    const polyMesh& mesh,
    const potential& pot,
    bool readFields
)
:
    Cloud<MoleculeType>(mesh, cloudName, false),
    moleculeCloud(),
    mesh_(mesh),
    pot_(pot),
    cellOccupancy_(mesh_.nCells()),
    il_(mesh_, pot_.pairPotentials().rCutMax(), false),
    constPropList_(),
    rndGen_(clock::getTime())
{
    if (readFields)
    {
        MoleculeType::readFields(*this);
    }

    buildConstProps();

    setSiteSizesAndPositions();

    removeHighEnergyOverlaps();

    calculateForce();
}


template<class MoleculeType>
Foam::MoleculeCloud<MoleculeType>::MoleculeCloud
(
    const word& cloudName,
    const polyMesh& mesh,
    const potential& pot,
    const IOdictionary& mdInitialiseDict,
    bool readFields
)
:
    Cloud<MoleculeType>(mesh, cloudName, false),
    moleculeCloud(),
    mesh_(mesh),
    pot_(pot),
    il_(mesh_, 0.0, false),
    constPropList_(),
    rndGen_(clock::getTime())
{
    if (readFields)
    {
        MoleculeType::readFields(*this);
    }

    this->clear();

    buildConstProps();

    initialiseMolecules(mdInitialiseDict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class MoleculeType>
void Foam::MoleculeCloud<MoleculeType>::evolve()
{
    typename MoleculeType::trackingData td0(*this, 0);
    Cloud<MoleculeType>::move(td0, mesh_.time().deltaTValue());

    typename MoleculeType::trackingData td1(*this, 1);
    Cloud<MoleculeType>::move(td1, mesh_.time().deltaTValue());

    typename MoleculeType::trackingData td2(*this, 2);
    Cloud<MoleculeType>::move(td2, mesh_.time().deltaTValue());

    calculateForce();

    typename MoleculeType::trackingData td3(*this, 3);
    Cloud<MoleculeType>::move(td3, mesh_.time().deltaTValue());

    info();
}


template<class MoleculeType>
void Foam::MoleculeCloud<MoleculeType>::calculateForce()
{
    buildCellOccupancy();

    // Set accumulated quantities to zero
    forAllIter(typename MoleculeCloud<MoleculeType>, *this, mol)
    {
        mol().siteForces() = vector::zero;

        mol().potentialEnergy() = 0.0;

        mol().rf() = tensor::zero;
    }

    calculatePairForce();

    calculateTetherForce();

    calculateExternalForce();
}


template<class MoleculeType>
void Foam::MoleculeCloud<MoleculeType>::info() const
{
    // Calculates and prints the mean momentum and energy in the system
    // and the number of molecules.

    vector totalLinearMomentum(vector::zero);

    vector totalAngularMomentum(vector::zero);

    scalar maxVelocityMag = 0.0;

    scalar totalMass = 0.0;

    scalar totalLinearKE = 0.0;

    scalar totalAngularKE = 0.0;

    scalar totalPE = 0.0;

    scalar totalrDotf = 0.0;

    //vector CentreOfMass(vector::zero);

    label nMols = this->size();

    label dofs = 0;

    {
        forAllConstIter(typename MoleculeCloud<MoleculeType>, *this, mol)
        {
            const label molId = mol().id();

            scalar molMass(this->constProps(molId).mass());

            totalMass += molMass;

            //CentreOfMass += mol().position()*molMass;
        }

        // if (nMols)
        // {
        //     CentreOfMass /= totalMass;
        // }

        forAllConstIter(typename MoleculeCloud<MoleculeType>, *this, mol)
        {
            const label molId = mol().id();

            const typename MoleculeType::constantProperties cP
            (
                this->constProps(molId)
            );

            scalar molMass(cP.mass());

            const diagTensor& molMoI(cP.momentOfInertia());

            const vector& molV(mol().v());

            const vector& molOmega(inv(molMoI) & mol().pi());

            vector molPiGlobal = mol().Q() & mol().pi();

            totalLinearMomentum += molV * molMass;

            totalAngularMomentum += molPiGlobal;
            //+((mol().position() - CentreOfMass) ^ (molV * molMass));

            if (mag(molV) > maxVelocityMag)
            {
                maxVelocityMag = mag(molV);
            }

            totalLinearKE += 0.5*molMass*magSqr(molV);

            totalAngularKE += 0.5*(molOmega & molMoI & molOmega);

            totalPE += mol().potentialEnergy();

            totalrDotf += tr(mol().rf());

            dofs += cP.degreesOfFreedom();
        }
    }

    scalar meshVolume = sum(mesh_.cellVolumes());

    if (Pstream::parRun())
    {
        reduce(totalLinearMomentum, sumOp<vector>());
        reduce(totalAngularMomentum, sumOp<vector>());
        reduce(maxVelocityMag, maxOp<scalar>());
        reduce(totalMass, sumOp<scalar>());
        reduce(totalLinearKE, sumOp<scalar>());
        reduce(totalAngularKE, sumOp<scalar>());
        reduce(totalPE, sumOp<scalar>());
        reduce(totalrDotf, sumOp<scalar>());
        reduce(nMols, sumOp<label>());
        reduce(dofs, sumOp<label>());
        reduce(meshVolume, sumOp<scalar>());
    }

    if (nMols)
    {
        Info<< nl << "Number of molecules in " << this->name() << " = "
            << nMols << nl
            << "    Overall number density = "
            << nMols/meshVolume << nl
            << "    Overall mass density = "
            << totalMass/meshVolume << nl
            << "    Average linear momentum per molecule = "
            << totalLinearMomentum/nMols << ' '
            << mag(totalLinearMomentum)/nMols << nl
            << "    Average angular momentum per molecule = "
            << totalAngularMomentum << ' '
            << mag(totalAngularMomentum)/nMols << nl
            << "    maximum |velocity| = "
            << maxVelocityMag << nl
            << "    Average linear KE per molecule = "
            << totalLinearKE/nMols << nl
            << "    Average angular KE per molecule = "
            << totalAngularKE/nMols << nl
            << "    Average PE per molecule = "
            << totalPE/nMols << nl
            << "    Average TE per molecule = "
            <<
            (
                  totalLinearKE
                + totalAngularKE
                + totalPE
            )
            /nMols
            << nl << endl;
    }
    else
    {
        Info<< nl << "No molecules in " << this->name() << nl << endl;
    }
}


template<class MoleculeType>
void Foam::MoleculeCloud<MoleculeType>::writeXYZ(const fileName& fName) const
{
    OFstream os(fName);

    os  << nSites() << nl
        << "MoleculeCloud<MoleculeType> site positions in angstroms" << nl;

    forAllConstIter(typename MoleculeCloud<MoleculeType>, *this, mol)
    {
        const typename MoleculeType::constantProperties& cP
        (
            constProps(mol().id())
        );

        forAll(mol().sitePositions(), i)
        {
            const point& sP = mol().sitePositions()[i];

            os  << pot_.siteIdList()[cP.sites()[i].siteId()]
                << ' ' << sP.x()*1e10
                << ' ' << sP.y()*1e10
                << ' ' << sP.z()*1e10
                << nl;
        }
    }
}


// ************************************************************************* //
