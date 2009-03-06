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
    Foam::ReactingMultiphaseParcel

Description
    Multiphase variant of the reacting parcel class with one/two-way coupling
    with the continuous phase.

SourceFiles
    ReactingMultiphaseParcelI.H
    ReactingMultiphaseParcel.C
    ReactingMultiphaseParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef ReactingMultiphaseParcel_H
#define ReactingMultiphaseParcel_H

#include "ReactingParcel.H"
#include "ReactingMultiphaseCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ParcelType>
class ReactingMultiphaseParcel;

template<class ParcelType>
Ostream& operator<<
(
    Ostream&,
    const ReactingMultiphaseParcel<ParcelType>&
);

/*---------------------------------------------------------------------------*\
                 Class ReactingMultiphaseParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class ReactingMultiphaseParcel
:
    public ReactingParcel<ParcelType>
{

public:

    //- Class to hold reacting particle constant properties
    class constantProperties
    :
        public ReactingParcel<ParcelType>::constantProperties
    {
        // Private data

            //- Latent heat of devolatilisation [J/kg]
            const scalar Ldevol_;

            //- Fraction of enthalpy retained by parcel due to surface
            //  reactions
            const scalar fh_;


    public:

        //- Constructor
        constantProperties(const dictionary& dict);

        // Access

            //- Return const access to the latent heat of devolatilisation
            inline scalar Ldevol() const;

            //- Return const access to the fraction of enthalpy retained by
            // parcel due to surface reactions
            inline scalar fh() const;
    };


    //- Class used to pass reacting tracking data to the trackToFace function
    class trackData
    :
        public ReactingParcel<ParcelType>::trackData
    {
        // Private data

            //- Reference to the cloud containing this particle
            ReactingMultiphaseCloud<ParcelType>& cloud_;

            //- Particle constant properties
            const constantProperties& constProps_;


    public:

        // Constructors

            //- Construct from components
            inline trackData
            (
                ReactingMultiphaseCloud<ParcelType>& cloud,
                const constantProperties& constProps,
                const interpolation<scalar>& rhoInterp,
                const interpolation<vector>& UInterp,
                const interpolation<scalar>& muInterp,
                const interpolation<scalar>& TInterp,
                const interpolation<scalar>& CpInterp,
                const interpolation<scalar>& pInterp,
                const vector& g
            );


        // Member functions

            //- Return access to the owner cloud
            inline ReactingMultiphaseCloud<ParcelType>& cloud();

            //- Return const access to the constant properties
            inline const constantProperties& constProps() const;
    };


private:

    // Private member functions

        //- Return the mixture effective specific heat capacity
        template<class TrackData>
        scalar cpEff(TrackData& td, const scalar p, const scalar T) const;

        //- Return the mixture effective enthalpy
        template<class TrackData>
        scalar HEff(TrackData& td, const scalar p, const scalar T) const;


protected:

    // Protected data

        // Parcel properties

            //- Mass fractions of gases []
            scalarField YGas_;

            //- Mass fractions of liquids []
            scalarField YLiquid_;

            //- Mass fractions of solids []
            scalarField YSolid_;

            //- Flag to say that the particle is allowed to combust
            //  Only true after volatile content falls below threshold value
            bool canCombust_;


    // Protected member functions

        //- Calculate Devolatilisation
        template<class TrackData>
        scalar calcDevolatilisation
        (
            TrackData& td,
            const scalar dt,
            const scalar T,
            scalarList& dMassMT
        );

        //- Calculate surface reactions
        template<class TrackData>
        void calcSurfaceReactions
        (
            TrackData& td,
            const scalar dt,
            const label cellI,
            const scalar T0,
            const scalar T1,
            const scalarList& dMassMT,
            scalarList& dMassSR,
            scalar& dhRet
        );


public:

    //- Runtime type information
    TypeName("ReactingMultiphaseParcel");

    friend class Cloud<ParcelType>;


    // Constructors

        //- Construct from components
        inline ReactingMultiphaseParcel
        (
            ReactingMultiphaseCloud<ParcelType>& owner,
            const label typeId,
            const vector& position,
            const label cellI,
            const scalar d0,
            const vector& U0,
            const scalar nParticle0,
            const scalarField& YGas0,
            const scalarField& YLiquid0,
            const scalarField& YSolid0,
            const scalarField& Y0,
            const constantProperties& constProps
        );

        //- Construct from Istream
        ReactingMultiphaseParcel
        (
            const Cloud<ParcelType>& c,
            Istream& is,
            bool readFields = true
        );

        //- Construct and return a clone
        autoPtr<ParcelType> clone() const
        {
            return autoPtr<ParcelType>
            (
                new ReactingMultiphaseParcel<ParcelType>(*this)
            );
        }


    // Member Functions

        // Access

            //- Return const access to mass fractions of gases
            inline const scalarField& YGas() const;

            //- Return const access to mass fractions of liquids
            inline const scalarField& YLiquid() const;

            //- Return const access to mass fractions of solids
            inline const scalarField& YSolid() const;


        // Edit

            //- Return access to mass fractions of gases
            inline scalarField& YGas();

            //- Return access to mass fractions of liquids
            inline scalarField& YLiquid();

            //- Return access to mass fractions of solids
            inline scalarField& YSolid();


        // Main calculation loop

            //- Update cell based quantities
            template<class TrackData>
            void updateCellQuantities
            (
                TrackData& td,
                const scalar dt,
                const label cellI
            );

            //- Update parcel properties over the time interval
            template<class TrackData>
            void calc
            (
                TrackData& td,
                const scalar dt,
                const label cellI
            );


        // I-O

            //- Read
            static void readFields(ReactingMultiphaseCloud<ParcelType>& c);

            //- Write
            static void writeFields
            (
                const ReactingMultiphaseCloud<ParcelType>& c
            );


    // Ostream Operator

        friend Ostream& operator<< <ParcelType>
        (
            Ostream&,
            const ReactingMultiphaseParcel<ParcelType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ReactingMultiphaseParcelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ReactingMultiphaseParcel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

