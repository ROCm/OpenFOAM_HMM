/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
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
    coalParcel

Description
    Coal parcel class

SourceFiles
    CoalParcel.C

\*---------------------------------------------------------------------------*/

#ifndef coalParcel_H
#define coalParcel_H

#include "ReactingMultiphaseParcel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class coalParcel Declaration
\*---------------------------------------------------------------------------*/

class coalParcel
:
    public ReactingMultiphaseParcel<coalParcel>
{
public:

    // Run-time type information
    TypeName("coalParcel");

    // Constructors

        //- Construct from owner, position, and cloud owner
        //  Other properties initialised as null
        coalParcel
        (
             ReactingMultiphaseCloud<coalParcel>& owner,
             const vector& position,
             const label cellI,
             const label tetFaceI,
             const label tetPtI
        );

        //- Construct from components
        coalParcel
        (
             ReactingMultiphaseCloud<coalParcel>& owner,
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
             const vector& pi0,
             const vector& tau0,
             const scalarField& YMixture0,
             const scalarField& YGas0,
             const scalarField& YLiquid0,
             const scalarField& YSolid0,
             const ReactingMultiphaseParcel<coalParcel>::constantProperties&
                 constProps
        );

        //- Construct from Istream
        coalParcel
        (
            const Cloud<coalParcel>& c,
            Istream& is,
            bool readFields = true
        );

        //- Construct as a copy
        coalParcel(const coalParcel& p);

        //- Construct as a copy
        coalParcel
        (
            const coalParcel& p,
            const ReactingMultiphaseCloud<coalParcel>& c
        );

        //- Construct and return a (basic particle) clone
        virtual autoPtr<Particle<coalParcel> > clone() const
        {
            return autoPtr<Particle<coalParcel> >(new coalParcel(*this));
        }

        //- Construct and return a (basic particle) clone
        virtual autoPtr<Particle<coalParcel> > clone
        (
            const Cloud<coalParcel>& c
        ) const
        {
            return autoPtr<Particle<coalParcel> >
            (
                new coalParcel
                (
                    *this,
                    static_cast<const ReactingMultiphaseCloud<coalParcel>&>(c)
                )
            );
        }


    //- Destructor
    virtual ~coalParcel();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
