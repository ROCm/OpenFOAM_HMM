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

Class
    Foam::PairCollision

Description

SourceFiles
    PairCollision.C

\*---------------------------------------------------------------------------*/

#ifndef PairCollision_H
#define PairCollision_H

#include "CollisionModel.H"
#include "InteractionLists.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class CloudType>
class PairModel;

template<class CloudType>
class WallModel;


/*---------------------------------------------------------------------------*\
                     Class PairCollision Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class PairCollision
:
    public CollisionModel<CloudType>
{
    // Static data

        //- Tolerance to determine flat wall interactions
        static scalar cosPhiMinFlatWall;

        //- Distance to multiply the distance to a wall with within
        //  which a duplicate flat wall interaction is defined.  If
        //  cosPhiFlatWall = 1 - e, this = sqrt(2e - e^2) < sqrt(3e)
        //  for all e < 1, giving a safe tolerance (i.e. no double
        //  interactions next to each other.)
        static scalar flatWallDuplicateExclusion;


    // Private data

        //- Cell occupancy information for each parcel
        List<DynamicList<typename CloudType::parcelType*> > cellOccupancy_;

        //- PairModel to calculate the interaction between two parcels
        autoPtr<PairModel<CloudType> > pairModel_;

        //- WallModel to calculate the interaction between the parcel and walls
        autoPtr<WallModel<CloudType> > wallModel_;

        //- Interactions lists determining which cells are in
        //  interaction range of each other
        InteractionLists<typename CloudType::parcelType> il_;


    // Private member functions

        //- Pre collision tasks
        void preInteraction();

        //- Interactions between real (on-processor) particles
        void realRealInteraction();

        //- Interactions between real and referred (off processor) particles
        void realReferredInteraction();

        //- Interactions with walls
        void wallInteraction();

        bool duplicatePointInList
        (
            const DynamicList<point>& existingPoints,
            const point& pointToTest,
            scalar duplicateRangeSqr
        ) const;

        bool duplicatePointInList
        (
            const DynamicList<point>& existingPoints,
            const point& pointToTest,
            const scalarList& duplicateRangeSqr
        ) const;

        //- Post collision tasks
        void postInteraction();

        //- Build the cell occupancy information for each parcel
        void buildCellOccupancy();

        //- Calculate the pair force between parcels
        void evaluatePair
        (
            typename CloudType::parcelType& pA,
            typename CloudType::parcelType& pB
        ) const;

        //- Calculate the wall forces on a parcel
        void evaluateWall
        (
            typename CloudType::parcelType& p,
            const List<point>& flatSites,
            const List<point>& sharpSites
        ) const;


public:

    //- Runtime type information
    TypeName("PairCollision");


    // Constructors

        //- Construct from dictionary
        PairCollision
        (
            const dictionary& dict,
            CloudType& owner
        );


    //- Destructor
    virtual ~PairCollision();


    // Member Functions

        //- Return the number of times to subcycle the current
        //  timestep to meet the criteria of the collision model.
        virtual label nSubCycles() const;

        //- Flag to indicate whether model activates injection model
        virtual bool active() const;

        // Collision function
        virtual void collide();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "PairCollision.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
