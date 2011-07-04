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

#include "polyatomicCloud.H"
#include "fvMesh.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTemplateTypeNameAndDebug(Cloud<polyatomic>, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::polyatomicCloud::buildConstProps()
{
    Info<< nl << "Reading mdProperties dictionary." << endl;

    const List<word>& idList(pot_.idList());

    constPropList_.setSize(idList.size());

    const List<word>& siteIdList(pot_.siteIdList());

    IOdictionary mdPropertiesDict
    (
        IOobject
        (
            "mdProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    );

    forAll(idList, i)
    {
        const word& id = idList[i];

        const dictionary& molDict(mdPropertiesDict.subDict(id));

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
                    "Foam::polyatomic::constantProperties::constantProperties"
                    "("
                    "const dictionary& dict"
                    ")"
                )
                    << siteId << " site not found."
                    << nl << abort(FatalError);
            }
        }

        constPropList_[i] = polyatomic::constantProperties
        (
            molDict,
            siteIds
        );
    }
}


void Foam::polyatomicCloud::setSiteSizesAndPositions()
{
    forAllIter(polyatomicCloud, *this, mol)
    {
        const polyatomic::constantProperties& cP = constProps(mol().id());

        mol().setSiteSizes(cP.nSites());

        mol().setSitePositions(cP);
    }
}


void Foam::polyatomicCloud::buildCellOccupancy()
{
    forAll(cellOccupancy_, cO)
    {
        cellOccupancy_[cO].clear();
    }

    forAllIter(polyatomicCloud, *this, mol)
    {
        cellOccupancy_[mol().cell()].append(&mol());
    }

    forAll(cellOccupancy_, cO)
    {
        cellOccupancy_[cO].shrink();
    }
}


void Foam::polyatomicCloud::calculatePairForce()
{
    PstreamBuffers pBufs(Pstream::nonBlocking);

    // Start sending referred data
    il_.sendReferredData(cellOccupancy(), pBufs);

    polyatomic* molI = NULL;
    polyatomic* molJ = NULL;

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
                    List<polyatomic*> cellJ =
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

        List<IDLList<polyatomic> >& referredMols = il_.referredParticles();

        forAll(ril, r)
        {
            const List<label>& realCells = ril[r];

            IDLList<polyatomic>& refMols = referredMols[r];

            forAllIter
            (
                IDLList<polyatomic>,
                refMols,
                refMol
            )
            {
                forAll(realCells, rC)
                {
                    List<polyatomic*> cellI = cellOccupancy_[realCells[rC]];

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


void Foam::polyatomicCloud::calculateTetherForce()
{
    const tetherPotentialList& tetherPot(pot_.tetherPotentials());

    forAllIter(polyatomicCloud, *this, mol)
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


void Foam::polyatomicCloud::calculateExternalForce()
{
    forAllIter(polyatomicCloud, *this, mol)
    {
        mol().a() += pot_.gravity();
    }
}


void Foam::polyatomicCloud::removeHighEnergyOverlaps()
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

    polyatomic* molI = NULL;
    polyatomic* molJ = NULL;

    {
        DynamicList<polyatomic*> molsToDelete;

        const labelListList& dil(il_.dil());

        forAll(dil, d)
        {
            forAll(cellOccupancy_[d], cellIMols)
            {
                molI = cellOccupancy_[d][cellIMols];

                forAll(dil[d], interactingCells)
                {
                    List<polyatomic*> cellJ =
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
        DynamicList<polyatomic*> molsToDelete;

        const labelListList& ril(il_.ril());

        List<IDLList<polyatomic> >& referredMols = il_.referredParticles();

        forAll(ril, r)
        {
            IDLList<polyatomic>& refMols = referredMols[r];

            forAllIter
            (
                IDLList<polyatomic>,
                refMols,
                refMol
            )
            {
                molJ = &refMol();

                const List<label>& realCells = ril[r];

                forAll(realCells, rC)
                {
                    label cellI = realCells[rC];

                    List<polyatomic*> cellIMols = cellOccupancy_[cellI];

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
                                // Remove one of the polyatomics
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

    Info<< tab << molsRemoved << " polyatomics removed" << endl;
}


void Foam::polyatomicCloud::initialiseMolecules
(
    const IOdictionary& mdInitialiseDict
)
{
    Info<< nl
        << "Initialising polyatomics in each zone specified in "
        << mdInitialiseDict.name()
        << endl;

    const cellZoneMesh& cellZones = mesh_.cellZones();

    if (!cellZones.size())
    {
        FatalErrorIn("void Foam::polyatomicCloud::initialiseMolecules")
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
                    FatalErrorIn("Foam::polyatomicCloud::initialiseMolecules")
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
                        WarningIn("polyatomicCloud::initialiseMolecules")
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

                        const polyatomic::constantProperties& cP
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
                        WarningIn("polyatomicCloud::initialiseMolecules")
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
                    FatalErrorIn("Foam::polyatomicCloud::initialiseMolecules")
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

                // Continue trying to place polyatomics as long as at
                // least one polyatomic is placed in each iteration.
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

                    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    // start of placement of polyatomics
                    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

                    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                    // end of placement of polyatomics
                    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                    if
                    (
                        totalZoneMols == 0
                     && !partOfLayerInBounds
                    )
                    {
                        WarningIn
                        (
                            "Foam::polyatomicCloud::initialiseMolecules()"
                        )
                            << "A whole layer of unit cells was placed "
                            << "outside the bounds of the mesh, but no "
                            << "polyatomics have been placed in zone '"
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


void Foam::polyatomicCloud::createMolecule
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
        FatalErrorIn("Foam::polyatomicCloud::createMolecule")
            << "Position specified does not correspond to a mesh cell." << nl
            << abort(FatalError);
    }

    point specialPosition(vector::zero);

    label special = 0;

    if (tethered)
    {
        specialPosition = position;

        special = polyatomic::SPECIAL_TETHERED;
    }

    const polyatomic::constantProperties& cP(constProps(id));

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
        new polyatomic
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


Foam::label Foam::polyatomicCloud::nSites() const
{
    label n = 0;

    forAllConstIter(polyatomicCloud, *this, mol)
    {
        n += constProps(mol().id()).nSites();
    }

    return n;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polyatomicCloud::polyatomicCloud
(
    const polyMesh& mesh,
    const potential& pot,
    bool readFields
)
:
    Cloud<polyatomic>(mesh, "polyatomicCloud", false),
    mesh_(mesh),
    pot_(pot),
    cellOccupancy_(mesh_.nCells()),
    il_(mesh_, pot_.pairPotentials().rCutMax(), false),
    constPropList_(),
    rndGen_(clock::getTime())
{
    if (readFields)
    {
        polyatomic::readFields(*this);
    }

    buildConstProps();

    setSiteSizesAndPositions();

    removeHighEnergyOverlaps();

    calculateForce();
}


Foam::polyatomicCloud::polyatomicCloud
(
    const polyMesh& mesh,
    const potential& pot,
    const IOdictionary& mdInitialiseDict,
    bool readFields
)
:
    Cloud<polyatomic>(mesh, "polyatomicCloud", false),
    mesh_(mesh),
    pot_(pot),
    il_(mesh_, 0.0, false),
    constPropList_(),
    rndGen_(clock::getTime())
{
    if (readFields)
    {
        polyatomic::readFields(*this);
    }

    clear();

    buildConstProps();

    initialiseMolecules(mdInitialiseDict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::polyatomicCloud::evolve()
{
    polyatomic::trackingData td0(*this, 0);
    Cloud<polyatomic>::move(td0, mesh_.time().deltaTValue());

    polyatomic::trackingData td1(*this, 1);
    Cloud<polyatomic>::move(td1, mesh_.time().deltaTValue());

    polyatomic::trackingData td2(*this, 2);
    Cloud<polyatomic>::move(td2, mesh_.time().deltaTValue());

    calculateForce();

    polyatomic::trackingData td3(*this, 3);
    Cloud<polyatomic>::move(td3, mesh_.time().deltaTValue());

    info();
}


void Foam::polyatomicCloud::calculateForce()
{
    buildCellOccupancy();

    // Set accumulated quantities to zero
    forAllIter(polyatomicCloud, *this, mol)
    {
        mol().siteForces() = vector::zero;

        mol().potentialEnergy() = 0.0;

        mol().rf() = tensor::zero;
    }

    calculatePairForce();

    calculateTetherForce();

    calculateExternalForce();
}


void Foam::polyatomicCloud::applyConstraintsAndThermostats
(
    const scalar targetTemperature,
    const scalar measuredTemperature
)
{
    scalar temperatureCorrectionFactor =
        sqrt(targetTemperature/max(VSMALL, measuredTemperature));

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

    forAllIter(polyatomicCloud, *this, mol)
    {
        mol().v() *= temperatureCorrectionFactor;

        mol().pi() *= temperatureCorrectionFactor;
    }
}


void Foam::polyatomicCloud::info() const
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
        forAllConstIter(polyatomicCloud, *this, mol)
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

        forAllConstIter(polyatomicCloud, *this, mol)
        {
            const label molId = mol().id();

            const polyatomic::constantProperties cP(this->constProps(molId));

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
        Info<< "Number of molecules in " << this->name() << " = "
            << nMols << nl
            << "Overall number density = "
            << nMols/meshVolume << nl
            << "Overall mass density = "
            << totalMass/meshVolume << nl
            << "Average linear momentum per molecule = "
            << totalLinearMomentum/nMols << ' '
            << mag(totalLinearMomentum)/nMols << nl
            << "Average angular momentum per molecule = "
            << totalAngularMomentum << ' '
            << mag(totalAngularMomentum)/nMols << nl
            << "maximum |velocity| = "
            << maxVelocityMag << nl
            << "Average linear KE per molecule = "
            << totalLinearKE/nMols << nl
            << "Average angular KE per molecule = "
            << totalAngularKE/nMols << nl
            << "Average PE per molecule = "
            << totalPE/nMols << nl
            << "Average TE per molecule = "
            <<
            (
                  totalLinearKE
                + totalAngularKE
                + totalPE
            )
            /nMols
            << endl;
    }
    else
    {
        Info<< "No molecules in " << this->name() << endl;
    }
}


void Foam::polyatomicCloud::writeXYZ(const fileName& fName) const
{
    OFstream os(fName);

    os  << nSites() << nl
        << "polyatomicCloud site positions in angstroms" << nl;

    forAllConstIter(polyatomicCloud, *this, mol)
    {
        const polyatomic::constantProperties& cP = constProps(mol().id());

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
