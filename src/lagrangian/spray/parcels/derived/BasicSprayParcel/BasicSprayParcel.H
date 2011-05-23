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
    Foam::BasicSprayParcel

Description

SourceFiles
    BasicSprayParcelI.H
    BasicSprayParcel.C
    BasicSprayParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef BasicSprayParcel_H
#define BasicSprayParcel_H

#include "SprayParcel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
template<class ThermoType>
class BasicSprayParcel;

/*---------------------------------------------------------------------------*\
                          Class BasicSprayParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class BasicSprayParcel
:
    public SprayParcel<BasicSprayParcel<ThermoType> >
{

public:

    //- The type of thermodynamics this parcel was instantiated for
    typedef ThermoType thermoType;

    //- Run-time type information
    TypeName("BasicSprayParcel");

    // Constructors

        //- Construct from owner, position, and cloud owner
        //  Other properties initialised as null
        BasicSprayParcel
        (
            SprayCloud<BasicSprayParcel>& owner,
            const vector& position,
            const label cellI
        );

        //- Construct from components
        BasicSprayParcel
        (
            SprayCloud<BasicSprayParcel>& owner,
            const vector& position,
            const label cellI,
            const label typeId,
            const scalar nParticle0,
            const scalar d0,
            const vector& U0,
            const scalarField& Y0,
            const scalar& dInit,
            const scalar& liquidCore,
            const scalar& KHindex,
            const scalar& y,
            const scalar& yDot,
            const scalar& tc,
            const scalar& ms,
            const scalar& injector,
            const scalar& tMom,
            const scalar& user,
            const typename SprayParcel<BasicSprayParcel>::
                constantProperties& constProps
        );

        //- Construct from Istream
        BasicSprayParcel
        (
            const Cloud<BasicSprayParcel>& c,
            Istream& is,
            bool readFields = true
        );

        //- Construct as a copy
        BasicSprayParcel(const BasicSprayParcel& p);

        //- Construct and return a clone
        autoPtr<BasicSprayParcel> clone() const
        {
            return
                autoPtr<BasicSprayParcel>
                (
                    new BasicSprayParcel(*this)
                );
        }


    //- Destructor
    virtual ~BasicSprayParcel();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "BasicSprayParcel.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
