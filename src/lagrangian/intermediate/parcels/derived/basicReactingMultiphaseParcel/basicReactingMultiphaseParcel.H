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
    Foam::BasicReactingMultiphaseParcel

Description


SourceFiles
    basicReactingMultiphaseParcel.C
    basicReactingMultiphaseParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef basicReactingMultiphaseParcel_H
#define basicReactingMultiphaseParcel_H

#include "ReactingMultiphaseParcel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class basicReactingMultiphaseParcel Declaration
\*---------------------------------------------------------------------------*/

class basicReactingMultiphaseParcel
:
    public ReactingMultiphaseParcel<basicReactingMultiphaseParcel>
{

public:

    //- Run-time type information
    TypeName("basicReactingMultiphaseParcel");

    // Constructors

        //- Construct from owner, position, and cloud owner
        //  Other properties initialised as null
        basicReactingMultiphaseParcel
        (
             ReactingMultiphaseCloud<basicReactingMultiphaseParcel>& owner,
             const vector& position,
             const label cellI,
             const label tetFaceI,
             const label tetPtI
        );

        //- Construct from components
        basicReactingMultiphaseParcel
        (
             ReactingMultiphaseCloud<basicReactingMultiphaseParcel>& owner,
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
             const scalarField& YGas0,
             const scalarField& YLiquid0,
             const scalarField& YSolid0,
             const scalarField& Y0,
             const ReactingMultiphaseParcel<basicReactingMultiphaseParcel>::
                constantProperties& constProps
        );

        //- Construct from Istream
        basicReactingMultiphaseParcel
        (
            const Cloud<basicReactingMultiphaseParcel>& c,
            Istream& is,
            bool readFields = true
        );

        //- Construct as a copy
        basicReactingMultiphaseParcel(const basicReactingMultiphaseParcel& p);

        //- Construct as a copy
        basicReactingMultiphaseParcel
        (
            const basicReactingMultiphaseParcel& p,
            const ReactingMultiphaseCloud<basicReactingMultiphaseParcel>& c
        );

        //- Construct and return a clone
        virtual autoPtr<Particle<basicReactingMultiphaseParcel> > clone() const
        {
            return autoPtr<Particle<basicReactingMultiphaseParcel> >
            (
                new basicReactingMultiphaseParcel(*this)
            );
        }

        //- Construct and return a clone
        virtual autoPtr<Particle<basicReactingMultiphaseParcel> > clone
        (
            const Cloud<basicReactingMultiphaseParcel>& c
        ) const
        {
            return autoPtr<Particle<basicReactingMultiphaseParcel> >
            (
                new basicReactingMultiphaseParcel
                (
                    *this,
                    static_cast
                    <
                        const ReactingMultiphaseCloud
                        <
                            basicReactingMultiphaseParcel
                        >&
                    >(c)
                )
            );
        }


    //- Destructor
    virtual ~basicReactingMultiphaseParcel();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
