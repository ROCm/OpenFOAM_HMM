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
    Foam::TrackedReactingParcel

Description
    Adds tracking to ReactingParcel

SourceFiles
    TrackedReactingParcelI.H
    TrackedReactingParcel.C
    TrackedReactingParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef TrackedReactingParcel_H
#define TrackedReactingParcel_H

#include "ReactingParcel.H"
#include "ParticleTrackingData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ParcelType>
class TrackedReactingParcel;

// Forward declaration of friend functions

template<class ParcelType>
Ostream& operator<<
(
    Ostream&,
    const TrackedReactingParcel<ParcelType>&
);

/*---------------------------------------------------------------------------*\
                  Class TrackedReactingParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class TrackedReactingParcel
:
    public ReactingParcel<ParcelType>,
    public ParticleTrackingData<ParcelType>
{
public:

    typedef typename ReactingParcel<ParcelType>::constantProperties
        constantProperties;

    // Static data

        //- Runtime type information
        TypeName("TrackedReactingParcel");


    // Constructors

        //- Construct from owner, position, and cloud owner
        //  Other properties initialised as null
        inline TrackedReactingParcel
        (
            ReactingCloud<ParcelType>& owner,
            const vector& position,
            const label cellI
        );

        //- Construct from components
        inline TrackedReactingParcel
        (
            ReactingCloud<ParcelType>& owner,
            const vector& position,
            const label cellI,
            const label typeId,
            const scalar nParticle0,
            const scalar d0,
            const vector& U0,
            const scalarField& Y0,
            const constantProperties& constProps
        );

        //- Construct from Istream
        TrackedReactingParcel
        (
            const Cloud<ParcelType>& c,
            Istream& is,
            bool readFields = true
        );

        //- Construct as a copy
        TrackedReactingParcel(const TrackedReactingParcel& p);

        //- Construct and return a clone
        autoPtr<TrackedReactingParcel> clone() const
        {
            return
                autoPtr<TrackedReactingParcel>
                (
                    new TrackedReactingParcel(*this)
                );
        }


    // Member Functions

        // I-O

            //- Read
            static void readFields(ReactingCloud<ParcelType>& c);

            //- Write
            static void writeFields(const ReactingCloud<ParcelType>& c);


    // Ostream Operator

        friend Ostream& operator<< <ParcelType>
        (
            Ostream&,
            const TrackedReactingParcel<ParcelType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "TrackedReactingParcelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "TrackedReactingParcel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

