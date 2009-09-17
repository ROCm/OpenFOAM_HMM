/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
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
    Foam::SpringSliderDashpot

Description
    Pair forces between particles colliding with a spring, slider, damper model

\*---------------------------------------------------------------------------*/

#ifndef SpringSliderDashpot_H
#define SpringSliderDashpot_H

#include "PairFunction.H"
#include "CollisionRecord.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                      Class SpringSliderDashpot Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class SpringSliderDashpot
:
    public PairFunction<CloudType>
{
    // Private data

        //- Effective Young's modulus value, assuming both particles have
        //  the same E value
        scalar Estar_;

        //- Effective shear modulus value, assuming both particles have
        //  the same Poisson's ratio and Young's modulus
        scalar Gstar_;

        //- Poisson's ratio of both particles
        scalar sigma_;

        //- alpha-coefficient, related to coefficient of restitution
        scalar alpha_;

        //- Spring power (b = 1 for linear, b = 3/2 for Hertzian)
        scalar b_;

        //- Coefficient of friction in for tangential sliding
        scalar mu_;

        //- The number of steps over which to resolve the minimum
        //  harmonic approximation of the collision period
        scalar collisionResolutionSteps_;


    // Private Member Functions

        //- Find the appropriate properties for determining the minimum
        //- allowable timestep
        void findMinMaxProperties
        (
            scalar& RMin,
            scalar& rhoMax,
            scalar& vMagMax
        ) const;

public:

    //- Runtime type information
    TypeName("SpringSliderDashpot");


    // Constructors

        //- Construct from dictionary
        SpringSliderDashpot(const dictionary& dict, CloudType& cloud);


    //- Destructor
    virtual ~SpringSliderDashpot();


    // Member Functions

        //- Whether the PairFunction has a timestep limit that will
        //  require subCycling
        virtual bool controlsTimestep() const;

        //- For PairFunctions that control the timestep, calculate the
        //  number of subCycles needed to satisfy the minimum
        //  allowable timestep
        virtual label nSubCycles() const;

        //- Calculate the pair interaction between parcels
        virtual void evaluatePair
        (
            typename CloudType::parcelType& pA,
            typename CloudType::parcelType& pB
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "SpringSliderDashpot.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
