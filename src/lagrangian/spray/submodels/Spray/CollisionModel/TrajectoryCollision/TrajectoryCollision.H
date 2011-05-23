/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2011 OpenCFD Ltd.
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

Class
    Foam::TrajectoryCollision

Description
    Trajectory collision model by N. Nordin.

\*---------------------------------------------------------------------------*/

#ifndef TrajectoryCollision_H
#define TrajectoryCollision_H

#include "CollisionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                       Class TrajectoryCollision Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class TrajectoryCollision
:
    public CollisionModel<CloudType>
{
private:
    // Reference to the coefficients dictionary
        dictionary coeffsDict_;

        scalar cSpace_;
        scalar cTime_;
        Switch coalescence_;
public:

    //- Runtime type information
    TypeName("Trajectory");


    // Constructors

        //- Construct from dictionary
        TrajectoryCollision
        (
            const dictionary& dict,
            CloudType& cloud
        );


    //- Destructor
    virtual ~TrajectoryCollision();


    // Member Functions

        virtual bool active() const;

        virtual bool update
        (
            const scalar& dt,
            Random& rndGen,
            vector& pos1,
            scalar& m1,
            scalar& d1,
            scalar& N1,
            vector& U,
            scalar& rho1,
            scalar& T1,
            scalarField& Y1,
            const scalar& sigma1,
            const label celli,
            const scalar voli,
            vector& pos2,
            scalar& m2,
            scalar& d2,
            scalar& N2,
            vector& U2,
            scalar& rho2,
            scalar& T2,
            scalarField& Y2,
            const scalar& sigma2,
            const label cellj,
            const scalar volj
        ) const;

        // 1 is the larger drop and 2 is the smaller
        bool collideSorted
        (
            const scalar& dt,
            Random& rndGen,
            vector& pos1,
            scalar& m1,
            scalar& d1,
            scalar& N1,
            vector& U,
            scalar& rho1,
            scalar& T1,
            scalarField& Y1,
            const scalar& sigma1,
            const label celli,
            const scalar voli,
            vector& pos2,
            scalar& m2,
            scalar& d2,
            scalar& N2,
            vector& U2,
            scalar& rho2,
            scalar& T2,
            scalarField& Y2,
            const scalar& sigma2,
            const label cellj,
            const scalar volj
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "TrajectoryCollision.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
