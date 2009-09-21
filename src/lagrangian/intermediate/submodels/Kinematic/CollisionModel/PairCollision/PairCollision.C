/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "PairCollision.H"
#include "PairModel.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
void Foam::PairCollision<CloudType>::preInteraction()
{
    // Set accumulated quantities to zero
    forAllIter(typename CloudType, this->owner(), iter)
    {
        typename CloudType::parcelType& p = iter();

        p.f() = vector::zero;

        p.tau() = vector::zero;
    }

    buildCellOccupancy();
}


template<class CloudType>
void Foam::PairCollision<CloudType>::realRealInteraction()
{
    const DirectInteractionList<typename CloudType::parcelType>& dil =
        il_.dil();

    typename CloudType::parcelType* pA_ptr = NULL;
    typename CloudType::parcelType* pB_ptr = NULL;

    forAll(dil, realCellI)
    {
        // Loop over all Parcels in cell A (a)
        forAll(cellOccupancy_[realCellI], a)
        {
            pA_ptr = cellOccupancy_[realCellI][a];

            forAll(dil[realCellI], interactingCells)
            {
                List<typename CloudType::parcelType*> cellBParcels =
                    cellOccupancy_[dil[realCellI][interactingCells]];

                // Loop over all Parcels in cell B (b)
                forAll(cellBParcels, b)
                {
                    pB_ptr = cellBParcels[b];

                    evaluatePair(*pA_ptr, *pB_ptr);
                }
            }

            // Loop over the other Parcels in cell A (aO)
            forAll(cellOccupancy_[realCellI], aO)
            {
                pB_ptr = cellOccupancy_[realCellI][aO];

                // Do not double-evaluate, compare pointers, arbitrary
                // order
                if (pB_ptr > pA_ptr)
                {
                    evaluatePair(*pA_ptr, *pB_ptr);
                }
            }
        }
    }
}


template<class CloudType>
void Foam::PairCollision<CloudType>::realReferredInteraction()
{
    ReferredCellList<typename CloudType::parcelType>& ril = il_.ril();

    // Loop over all referred cells
    forAll(ril, refCellI)
    {
        ReferredCell<typename CloudType::parcelType>& refCell = ril[refCellI];

        const labelList& realCells = refCell.realCellsForInteraction();

        // Loop over all referred parcels in the referred cell

        forAllIter
        (
            typename IDLList<typename CloudType::parcelType>,
            refCell,
            referredParcel
        )
        {
            // Loop over all real cells in that the referred cell is
            // to supply interactions to

            forAll(realCells, realCellI)
            {
                List<typename CloudType::parcelType*> realCellParcels =
                    cellOccupancy_[realCells[realCellI]];

                forAll(realCellParcels, realParcelI)
                {
                    evaluatePair
                    (
                        *realCellParcels[realParcelI],
                        referredParcel()
                    );
                }
            }
        }
    }
}


template<class CloudType>
void Foam::PairCollision<CloudType>::wallInteraction()
{
    const polyMesh& mesh = this->owner().mesh();

    const DirectInteractionList<typename CloudType::parcelType>& dil =
        il_.dil();

    ReferredCellList<typename CloudType::parcelType>& ril = il_.ril();

    DynamicList<point> allWallInteractionSites;
    DynamicList<point> flatWallInteractionSites;
    DynamicList<point> sharpWallInteractionSites;

    forAll(dil, realCellI)
    {
        // The real wall faces in range of this real cell
        const labelList& realWallFaces = dil.wallFaces()[realCellI];

        // The labels of referred cells in range of this real cell
        const labelList& referredCellsInRange =
        dil.referredCellsForInteraction()[realCellI];

        // Loop over all Parcels in cell
        forAll(cellOccupancy_[realCellI], cellParticleI)
        {
            allWallInteractionSites.clear();
            flatWallInteractionSites.clear();
            sharpWallInteractionSites.clear();

            typename CloudType::parcelType& p =
                *cellOccupancy_[realCellI][cellParticleI];

            const point& pt = p.position();

            // real wallFace interactions

            forAll(realWallFaces, realWallFaceI)
            {
                label realFaceI = realWallFaces[realWallFaceI];

                pointHit nearest = mesh.faces()[realFaceI].nearestPoint
                (
                    pt,
                    mesh.points()
                );

                if (nearest.distance() < p.r())
                {
                    vector normal = mesh.faceAreas()[realFaceI];

                    normal /= mag(normal);

                    vector pW = nearest.rawPoint() - pt;

                    scalar normalAlignment = normal & pW/mag(pW);

                    allWallInteractionSites.append(nearest.rawPoint());

                    if (normalAlignment > 1 - SMALL)
                    {
                        flatWallInteractionSites.append(nearest.rawPoint());
                    }
                }
            }

            // referred wallFace interactions

            forAll(referredCellsInRange, refCellInRangeI)
            {
                ReferredCell<typename CloudType::parcelType>& refCell =
                    ril[referredCellsInRange[refCellInRangeI]];

                const labelList& refWallFaces = refCell.wallFaces();

                forAll(refWallFaces, refWallFaceI)
                {
                    label refFaceI = refWallFaces[refWallFaceI];

                    pointHit nearest = refCell.faces()[refFaceI].nearestPoint
                    (
                        pt,
                        refCell.points()
                    );

                    if (nearest.distance() < p.r())
                    {
                        vector normal = refCell.faceAreas()[refFaceI];

                        normal /= mag(normal);

                        vector pW = nearest.rawPoint() - pt;

                        scalar normalAlignment = normal & pW/mag(pW);

                        allWallInteractionSites.append(nearest.rawPoint());

                        if (normalAlignment > 1 - SMALL)
                        {
                            flatWallInteractionSites.append(nearest.rawPoint());
                        }
                    }
                }
            }

            Pout<< flatWallInteractionSites << endl;

            forAll(flatWallInteractionSites, siteI)
            {

                scalar nu = this->owner().constProps().poissonsRatio();

                scalar E = this->owner().constProps().youngsModulus();

                scalar b = 1.5;

                scalar alpha = 0.52;

                vector r_PW = pt - flatWallInteractionSites[siteI];

                scalar normalOverlapMag = p.r() - mag(r_PW);

                vector rHat_PW = r_PW/(mag(r_PW) + VSMALL);

                scalar kN = (4.0/3.0)*sqrt(p.r())*E/(2.0*(1.0 - sqr(nu)));

                scalar etaN = alpha*sqrt(p.mass()*kN)*pow025(normalOverlapMag);

                vector fN_PW =
                    rHat_PW
                   *(kN*pow(normalOverlapMag, b) - etaN*(p.U() & rHat_PW));

                p.f() += fN_PW;

                Pout<< "Wall force " << fN_PW << endl;
            }
        }
    }
}


template<class CloudType>
void Foam::PairCollision<CloudType>::postInteraction()
{
    // Delete any collision records where no collision occurred this step

    Info<< "    Update collision records" << endl;

    forAllIter(typename CloudType, this->owner(), iter)
    {
        typename CloudType::parcelType& p = iter();

        p.collisionRecords().update();
    }
}


template<class CloudType>
void Foam::PairCollision<CloudType>::buildCellOccupancy()
{
    Info<< "    Build cell occupancy" << endl;

    forAll(cellOccupancy_, cO)
    {
        cellOccupancy_[cO].clear();
    }

    forAllIter(typename CloudType, this->owner(), iter)
    {
        cellOccupancy_[iter().cell()].append(&iter());
    }

    il_.ril().referParticles(cellOccupancy_);
}


template<class CloudType>
void Foam::PairCollision<CloudType>::evaluatePair
(
    typename CloudType::parcelType& pA,
    typename CloudType::parcelType& pB
) const
{
    pairModel_->evaluatePair(pA, pB);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PairCollision<CloudType>::PairCollision
(
    const dictionary& dict,
    CloudType& owner
)
:
    CollisionModel<CloudType>(dict, owner, typeName),
    cellOccupancy_(owner.mesh().nCells()),
    pairModel_
    (
        PairModel<CloudType>::New
        (
            this->coeffDict(),
            this->owner()
        )
    ),
    il_
    (
        owner.mesh(),
        sqr(readScalar(this->coeffDict().lookup("maxInteractionDistance"))),
        true
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PairCollision<CloudType>::~PairCollision()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::label Foam::PairCollision<CloudType>::nSubCycles() const
{
    if (pairModel_->controlsTimestep())
    {
        label nSubCycles = returnReduce
        (
            pairModel_->nSubCycles(), maxOp<label>()
        );

        if(nSubCycles > 1)
        {
            Info<< nSubCycles << " move-collide subCycles" << endl;
        }

        return nSubCycles;
    }
    else
    {
        return 1;
    }
}


template<class CloudType>
bool Foam::PairCollision<CloudType>::active() const
{
    return true;
}


template<class CloudType>
void Foam::PairCollision<CloudType>::collide()
{
    Info<< "Calculating collisions" << endl;

    preInteraction();

    realRealInteraction();

    realReferredInteraction();

    wallInteraction();

    postInteraction();
}


// ************************************************************************* //
