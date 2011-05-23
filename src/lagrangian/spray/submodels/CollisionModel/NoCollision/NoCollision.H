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
    Foam::NoCollision

Description
    Dummy collision model for 'none'

\*---------------------------------------------------------------------------*/

#ifndef NoCollision_H
#define NoCollision_H

#include "CollisionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                       Class NoCollision Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoCollision
:
    public CollisionModel<CloudType>
{
public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct from dictionary
        NoCollision(const dictionary& dict, CloudType& cloud);

        //- Construct copy
        NoCollision(const NoCollision<CloudType>& cm);

        //- Construct and return a clone
        virtual autoPtr<CollisionModel<CloudType> > clone() const
        {
            return autoPtr<CollisionModel<CloudType> >
            (
                new NoCollision<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~NoCollision();


    // Member Functions

        virtual bool update
        (
            const scalar dt,
            cachedRandom& rndGen,
            vector& pos1,
            scalar& m1,
            scalar& d1,
            scalar& N1,
            vector& U,
            scalar& rho1,
            scalar& T1,
            scalarField& Y1,
            const scalar sigma1,
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
            const scalar sigma2,
            const label cellj,
            const scalar volj
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "NoCollision.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
