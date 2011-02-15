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

Class
    Foam::basicKinematicParcel

Description


SourceFiles
    basicKinematicParcel.C

\*---------------------------------------------------------------------------*/

#ifndef basicKinematicParcel_H
#define basicKinematicParcel_H

#include "KinematicParcel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class basicKinematicParcel Declaration
\*---------------------------------------------------------------------------*/

class basicKinematicParcel
:
    public KinematicParcel<basicKinematicParcel>
{

public:

    //- Run-time type information
    TypeName("basicKinematicParcel");


    // Constructors

        //- Construct from owner, position, and cloud owner
        //  Other properties initialised as null
        basicKinematicParcel
        (
            KinematicCloud<basicKinematicParcel>& owner,
            const vector& position,
            const label cellI,
            const label tetFaceI,
            const label tetPtI
        );

        //- Construct from components
        basicKinematicParcel
        (
            KinematicCloud<basicKinematicParcel>& owner,
            const vector& position,
            const label cellI,
            const label tetFaceI,
            const label tetPtI,
            const label typeId,
            const scalar nParticle0,
            const scalar d0,
            const scalar dTarget0,
            const vector& U0,
            const vector& f0,
            const vector& angularMomentum0,
            const vector& torque0,
            const constantProperties& constProps
        );

        //- Construct from Istream
        basicKinematicParcel
        (
            const Cloud<basicKinematicParcel>& c,
            Istream& is,
            bool readFields = true
        );

        //- Construct as a copy
        basicKinematicParcel(const basicKinematicParcel& p);

        //- Construct as a copy
        basicKinematicParcel
        (
            const basicKinematicParcel& p,
            const KinematicCloud<basicKinematicParcel>& c
        );

        //- Construct and return a clone
        virtual autoPtr<Particle<basicKinematicParcel> > clone() const
        {
            return autoPtr<Particle<basicKinematicParcel> >
            (
                new basicKinematicParcel(*this)
            );
        }

        //- Construct and return a clone
        virtual autoPtr<Particle<basicKinematicParcel> > clone
        (
            const Cloud<basicKinematicParcel>& c
        ) const
        {
            return autoPtr<Particle<basicKinematicParcel> >
            (
                new basicKinematicParcel
                (
                    *this,
                    static_cast<const KinematicCloud<basicKinematicParcel>&>(c)
                )
            );
        }


    //- Destructor
    virtual ~basicKinematicParcel();
};


template<>
inline bool contiguous<basicKinematicParcel>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
