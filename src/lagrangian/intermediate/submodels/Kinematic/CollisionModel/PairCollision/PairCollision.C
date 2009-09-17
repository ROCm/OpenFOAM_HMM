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
#include "PairFunction.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

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
    pairFunction_->evaluatePair(pA, pB);
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
    pairFunction_
    (
        PairFunction<CloudType>::New
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
    if (pairFunction_->controlsTimestep())
    {
        label nSubCycles = returnReduce
        (
            pairFunction_->nSubCycles(), maxOp<label>()
        );

        Info<< nSubCycles << " move-collide subCycles" << endl;

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

    // Set accumulated quantities to zero
    forAllIter(typename CloudType, this->owner(), iter)
    {
        typename CloudType::parcelType& p = iter();

        p.f() = vector::zero;

        p.tau() = vector::zero;
    }

    buildCellOccupancy();

    const DirectInteractionList<typename CloudType::parcelType>& dil(il_.dil());

    typename CloudType::parcelType* pA_ptr = NULL;
    typename CloudType::parcelType* pB_ptr = NULL;

    // real-real interactions

    forAll(dil, d)
    {
        // Loop over all Parcels in cell A (a)
        forAll(cellOccupancy_[d], a)
        {
            pA_ptr = cellOccupancy_[d][a];

            forAll(dil[d], interactingCells)
            {
                List<typename CloudType::parcelType*> cellBParcels =
                    cellOccupancy_[dil[d][interactingCells]];

                // Loop over all Parcels in cell B (b)
                forAll(cellBParcels, b)
                {
                    pB_ptr = cellBParcels[b];

                    evaluatePair(*pA_ptr, *pB_ptr);
                }
            }

            // Loop over the other Parcels in cell A (aO)
            forAll(cellOccupancy_[d], aO)
            {
                pB_ptr = cellOccupancy_[d][aO];

                // Do not double-evaluate, compare pointers, arbitrary
                // order
                if (pB_ptr > pA_ptr)
                {
                    evaluatePair(*pA_ptr, *pB_ptr);
                }
            }
        }
    }

    // real-referred interactions

    ReferredCellList<typename CloudType::parcelType>& ril(il_.ril());

    // Loop over all referred cells
    forAll(ril, refCellI)
    {
        ReferredCell<typename CloudType::parcelType>& refCell =
            ril[refCellI];

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

    Info<< "    ADD COLLISIONS WITH WALLS HERE" << endl;
    //     << " DOES NOT NEED TO BE A TRACKING OPERATION."
    //     << " CALCULATE DISTANCE TO SURFACES OF WALL TYPE AND APPLY "
    //     << "WALL FORCE MODEL" << endl;

    // Delete any collision records where no collision occurred this step

    Info<< "    Update collision records" << endl;

    forAllIter(typename CloudType, this->owner(), iter)
    {
        typename CloudType::parcelType& p = iter();

        p.collisionRecords().update();
    }
}


// ************************************************************************* //
