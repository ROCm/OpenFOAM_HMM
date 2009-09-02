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
    Foam::DeterministicPairForce

Description

SourceFiles
    DeterministicPairForce.C

\*---------------------------------------------------------------------------*/

#ifndef DeterministicPairForce_H
#define DeterministicPairForce_H

#include "CollisionModel.H"
#include "InteractionLists.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class DeterministicPairForce Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class DeterministicPairForce
:
    public CollisionModel<CloudType>
{
    // Private data

        //- Cell occupancy information for each parcel
        List<DynamicList<typename CloudType::parcelType*> > cellOccupancy_;

        //- Interactions lists determining which cells are in
        //  interaction range of each other
        InteractionLists<typename CloudType::parcelType> il_;


    // Private member functions

        //- Build the cell occupancy information for each parcel
        void buildCellOccupancy();


protected:

    // Protected member functions

        //- Calculate the pair force between parcels
        virtual void evaluatePair
        (
            typename CloudType::parcelType& pA,
            typename CloudType::parcelType& pB
        ) const;


public:

    //- Runtime type information
    TypeName("DeterministicPairForce");


    // Constructors

        //- Construct from dictionary
        DeterministicPairForce
        (
            const dictionary& dict,
            CloudType& owner
        );


    //- Destructor
    virtual ~DeterministicPairForce();


    // Member Functions

        //- Flag to indicate whether model activates injection model
        virtual bool active() const;

        // Collision function
        virtual void collide();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "DeterministicPairForce.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
