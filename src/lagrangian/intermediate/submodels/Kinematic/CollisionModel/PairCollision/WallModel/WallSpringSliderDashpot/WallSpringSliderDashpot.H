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
    Foam::WallSpringSliderDashpot

Description
    Forces between particles and walls, interacting with a spring,
    slider, damper model

\*---------------------------------------------------------------------------*/

#ifndef WallSpringSliderDashpot_H
#define WallSpringSliderDashpot_H

#include "WallModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                    Class WallSpringSliderDashpot Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class WallSpringSliderDashpot
:
    public WallModel<CloudType>
{
    // Private data

        //- Young's modulus of the wall
        scalar E_;

        //- Poisson's ratio of the wall
        scalar nu_;

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
    TypeName("WallSpringSliderDashpot");


    // Constructors

        //- Construct from dictionary
        WallSpringSliderDashpot(const dictionary& dict, CloudType& cloud);


    //- Destructor
    virtual ~WallSpringSliderDashpot();


    // Member Functions

        //- Whether the WallModel has a timestep limit that will
        //  require subCycling
        virtual bool controlsTimestep() const;

        //- For WallModels that control the timestep, calculate the
        //  number of subCycles needed to satisfy the minimum
        //  allowable timestep
        virtual label nSubCycles() const;

        //- Calculate the wall interaction between parcels
        virtual void evaluateWall
        (
            typename CloudType::parcelType& p,
            const List<point>& flatSites,
            const List<point>& sharpSites
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "WallSpringSliderDashpot.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
