/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2010 OpenCFD Ltd.
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
    Foam::basicThermoParcel

Description


SourceFiles
    basicThermoParcel.C

\*---------------------------------------------------------------------------*/

#ifndef basicThermoParcel_H
#define basicThermoParcel_H

#include "ThermoParcel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class basicThermoParcel Declaration
\*---------------------------------------------------------------------------*/

class basicThermoParcel
:
    public ThermoParcel<basicThermoParcel>
{

public:

    //- Runtime type information
    TypeName("basicThermoParcel");

    // Constructors

       //- Construct from owner, position, and cloud owner
       //  Other properties initialised as null
       basicThermoParcel
       (
            ThermoCloud<basicThermoParcel>& owner,
            const vector position,
            const label cellI,
            const label tetFaceI,
            const label tetPtI
       );

       //- Construct from components
       basicThermoParcel
       (
            ThermoCloud<basicThermoParcel>& owner,
            const vector position,
            const label cellI,
            const label tetFaceI,
            const label tetPtI,
            const label typeId,
            const scalar nParticle0,
            const scalar d0,
            const scalar dTarget0,
            const vector U0,
            const vector& f0,
            const vector& angularMomentum0,
            const vector& torque0,
            const constantProperties& constProps
        );

        //- Construct from Istream
        basicThermoParcel
        (
            const Cloud<basicThermoParcel>& c,
            Istream& is,
            bool readFields = true
        );

        //- Construct as a copy
        basicThermoParcel(const basicThermoParcel& p);

        //- Construct as a copy
        basicThermoParcel
        (
            const basicThermoParcel& p,
            const ThermoCloud<basicThermoParcel>& c
        );

        //- Construct and return a clone
        virtual autoPtr<Particle<basicThermoParcel> > clone() const
        {
            return autoPtr<Particle<basicThermoParcel> >
            (
                new basicThermoParcel(*this)
            );
        }

        //- Construct and return a clone
        virtual autoPtr<Particle<basicThermoParcel> > clone
        (
            const Cloud<basicThermoParcel>& c
        ) const
        {
            return autoPtr<Particle<basicThermoParcel> >
            (
                new basicThermoParcel
                (
                    *this,
                    static_cast<const ThermoCloud<basicThermoParcel>&>(c)
                )
            );
        }


    //- Destructor
    virtual ~basicThermoParcel();
};


template<>
inline bool contiguous<basicThermoParcel>()
{
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
