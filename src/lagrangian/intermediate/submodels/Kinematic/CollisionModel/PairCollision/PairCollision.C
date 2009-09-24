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
#include "WallModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::PairCollision<CloudType>::cosPhiMinFlatWall = 1 - SMALL;

template<class CloudType>
Foam::scalar Foam::PairCollision<CloudType>::flatWallDuplicateExclusion =
    sqrt(3*SMALL);


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

    const ReferredCellList<typename CloudType::parcelType>& ril = il_.ril();

    // Storage for the wall interaction sites
    DynamicList<point> flatSites;
    DynamicList<scalar> flatSiteExclusionDistancesSqr;
    DynamicList<point> otherSites;
    DynamicList<scalar> otherSiteDistances;
    DynamicList<point> sharpSites;
    DynamicList<scalar> sharpSiteExclusionDistancesSqr;

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
            flatSites.clear();
            flatSiteExclusionDistancesSqr.clear();
            otherSites.clear();
            otherSiteDistances.clear();
            sharpSites.clear();
            sharpSiteExclusionDistancesSqr.clear();

            typename CloudType::parcelType& p =
                *cellOccupancy_[realCellI][cellParticleI];

            const point& pos = p.position();

            // real wallFace interactions

            forAll(realWallFaces, realWallFaceI)
            {
                label realFaceI = realWallFaces[realWallFaceI];

                pointHit nearest = mesh.faces()[realFaceI].nearestPoint
                (
                    pos,
                    mesh.points()
                );

                if (nearest.distance() < p.r())
                {
                    vector normal = mesh.faceAreas()[realFaceI];

                    normal /= mag(normal);

                    const vector& nearPt = nearest.rawPoint();

                    vector pW = nearPt - pos;

                    scalar normalAlignment = normal & pW/mag(pW);

                    if (normalAlignment > cosPhiMinFlatWall)
                    {
                        // Guard against a flat interaction being
                        // present on the boundary of two or more
                        // faces, which would create duplicate contact
                        // points. Duplicates are discarded.
                        if
                        (
                            !duplicatePointInList
                            (
                                flatSites,
                                nearPt,
                                sqr(p.r()*flatWallDuplicateExclusion)
                            )
                        )
                        {
                            flatSites.append(nearPt);

                            flatSiteExclusionDistancesSqr.append
                            (
                                sqr(p.r()) - sqr(nearest.distance())
                            );
                        }
                    }
                    else
                    {
                        otherSites.append(nearPt);

                        otherSiteDistances.append(nearest.distance());
                    }
                }
            }

            // referred wallFace interactions

            forAll(referredCellsInRange, refCellInRangeI)
            {
                const ReferredCell<typename CloudType::parcelType>& refCell =
                    ril[referredCellsInRange[refCellInRangeI]];

                const labelList& refWallFaces = refCell.wallFaces();

                forAll(refWallFaces, refWallFaceI)
                {
                    label refFaceI = refWallFaces[refWallFaceI];

                    pointHit nearest = refCell.faces()[refFaceI].nearestPoint
                    (
                        pos,
                        refCell.points()
                    );

                    if (nearest.distance() < p.r())
                    {
                        vector normal = refCell.faceAreas()[refFaceI];

                        normal /= mag(normal);

                        const vector& nearPt = nearest.rawPoint();

                        vector pW = nearPt - pos;

                        scalar normalAlignment = normal & pW/mag(pW);

                        if (normalAlignment > cosPhiMinFlatWall)
                        {
                            // Guard against a flat interaction being
                            // present on the boundary of two or more
                            // faces, which would create duplicate contact
                            // points. Duplicates are discarded.
                            if
                            (
                                !duplicatePointInList
                                (
                                    flatSites,
                                    nearPt,
                                    sqr(p.r()*flatWallDuplicateExclusion)
                                )
                            )
                            {
                                flatSites.append(nearPt);

                                flatSiteExclusionDistancesSqr.append
                                (
                                    sqr(p.r()) - sqr(nearest.distance())
                                );
                            }
                        }
                        else
                        {
                            otherSites.append(nearPt);

                            otherSiteDistances.append(nearest.distance());
                        }
                    }
                }
            }

            // All flat interaction sites found, now classify the
            // other sites as being in range of a flat interaction, or
            // a sharp interaction, being aware of not duplicating the
            // sharp interaction sites.

            // The "other" sites need to evaluated in order of
            // ascending distance to their nearest point so that
            // grouping occurs around the closest in any group

            labelList sortedOtherSiteIndices;

            sortedOrder(otherSiteDistances, sortedOtherSiteIndices);

            forAll(sortedOtherSiteIndices, siteI)
            {
                label orderedIndex = sortedOtherSiteIndices[siteI];

                const point& otherPt = otherSites[orderedIndex];

                if
                (
                    !duplicatePointInList
                    (
                        flatSites,
                        otherPt,
                        flatSiteExclusionDistancesSqr
                    )
                )
                {
                    // Not in range of a flat interaction, must be a
                    // sharp interaction.

                    if
                    (
                        !duplicatePointInList
                        (
                            sharpSites,
                            otherPt,
                            sharpSiteExclusionDistancesSqr
                        )
                    )
                    {
                        sharpSites.append(otherPt);

                        sharpSiteExclusionDistancesSqr.append
                        (
                            sqr(p.r()) - sqr(otherSiteDistances[orderedIndex])
                        );
                    }
                }
            }

            evaluateWall(p, flatSites, sharpSites);
        }
    }
}


template<class CloudType>
bool Foam::PairCollision<CloudType>::duplicatePointInList
(
    const DynamicList<point>& existingPoints,
    const point& pointToTest,
    scalar duplicateRangeSqr
) const
{
    forAll(existingPoints, i)
    {
        if (magSqr(existingPoints[i] - pointToTest) < duplicateRangeSqr)
        {
            return true;
        }
    }

    return false;
}


template<class CloudType>
bool Foam::PairCollision<CloudType>::duplicatePointInList
(
    const DynamicList<point>& existingPoints,
    const point& pointToTest,
    const scalarList& duplicateRangeSqr
) const
{
    forAll(existingPoints, i)
    {
        if (magSqr(existingPoints[i] - pointToTest) < duplicateRangeSqr[i])
        {
            return true;
        }
    }

    return false;
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


template<class CloudType>
void Foam::PairCollision<CloudType>::evaluateWall
(
    typename CloudType::parcelType& p,
    const List<point>& flatSites,
    const List<point>& sharpSites
) const
{
    wallModel_->evaluateWall(p, flatSites, sharpSites);
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
    wallModel_
    (
        WallModel<CloudType>::New
        (
            this->coeffDict(),
            this->owner()
        )
    ),
    il_
    (
        owner.mesh(),
        sqr(readScalar(this->coeffDict().lookup("maxInteractionDistance"))),
        false
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
