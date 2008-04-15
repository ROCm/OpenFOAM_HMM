/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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
    Foam::ReactingParcel

Description
    Reacting parcel class with one-way coupling with the continuous
    phase.

    Includes thermo parcel sub-models, plus:
    - combustion

SourceFiles
    reactingParcelI.H
    ReactingParcel.C
    reactingParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef ReactingParcel_H
#define ReactingParcel_H

#include "IOstream.H"
#include "autoPtr.H"
#include "interpolationCellPoint.H"
#include "contiguous.H"

#include "ThermoParcel.H"
#include "ReactingCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ParcelType>
class ReactingParcel;

template<class ParcelType>
Ostream& operator<<
(
    Ostream&,
    const ReactingParcel<ParcelType>&
);

/*---------------------------------------------------------------------------*\
                        Class ReactingParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class ReactingParcel
:
    public ThermoParcel<ParcelType>
{

    // Private data

        //- Initial particle mass
        scalar mass0_;

        //- Mass fractions of mixture
        scalarField YMixture_;

        //- Mass fractions of gases
        scalarField YGas_;

        //- Mass fractions of liquids
        scalarField YLiquid_;

        //- Mass fractions of solids
        scalarField YSolid_;

        //- Flag to say that the particle is allowed to combust
        //  Only true one ALL volatiles have been evolved from the particle
        bool canCombust_;


protected:

    // Protected member functions

        template<class TrackingData>
        void calcCoupled
        (
            TrackingData& td,
            const label celli,
            const scalar dt,
            const scalar rhoc,
            vector& Uc,
            const scalar muc,
            const scalar Tc,
            const scalar Cpc,
            const scalar pc
        );

        template<class TrackingData>
        void calcUncoupled
        (
            TrackingData& td,
            const label celli,
            const scalar dt,
            const scalar rhoc,
            vector& Uc,
            const scalar muc,
            const scalar Tc,
            const scalar Cpc,
            const scalar pc
        );

        //- Calculate mass transfer
        template<class TrackingData>
        void calcMassTransfer
        (
            TrackingData& td,
            const scalar dt,
            const scalar T0,
            const scalar T1,
            scalarList& dMassMT
        );

        //- Calculate surface reactions
        template<class TrackingData>
        void calcSurfaceReactions
        (
            TrackingData& td,
            const scalar dt,
            const label celli,
            const scalar rhoc,
            const scalar Tc,
            const scalar T0,
            const scalar T1,
            scalar& dMassMTSR,
            scalarList& dMassMT
        );


public:

    //- Runtime type information
    TypeName("ReactingParcel");

    friend class Cloud<ParcelType>;


    //- Class to hold particle constant properties
    class constantProperties
    :
        public ThermoParcel<ParcelType>::constantProperties
    {

        // Private data

            //- Vapourisation temperature
            const scalar Tvap_;

            //- Boiling point
            const scalar Tbp_;


    public:

        // Constructors
        constantProperties
        (
            const dictionary& dict
        );

        // Member functions

            // Access

            inline const scalar Tvap() const;
            inline const scalar Tbp() const;
    };


    //- Class used to pass tracking data to the trackToFace function
    class trackData
    :
        public ThermoParcel<ParcelType>::trackData
    {

        // Private data

            //- Reference to the cloud containing this particle
            ReactingCloud<ParcelType>& cloud_;

            //- Particle constant properties
            const constantProperties& constProps_;

            //- Interpolators for continuous phase fields

                const interpolationCellPoint<scalar>& pInterp_;

    public:

        // Constructors

            inline trackData
            (
                ReactingCloud<ParcelType>& cloud,
                const constantProperties& constProps,
                const interpolationCellPoint<scalar>& rhoInterp,
                const interpolationCellPoint<vector>& UInterp,
                const interpolationCellPoint<scalar>& muInterp,
                const interpolationCellPoint<scalar>& TInterp,
                const interpolationCellPoint<scalar>& CpInterp,
                const interpolationCellPoint<scalar>& pInterp,
                const vector& g
            );


        // Member functions

            inline ReactingCloud<ParcelType>& cloud();

            inline const constantProperties& constProps() const;

            inline const interpolationCellPoint<scalar>& pInterp() const;
    };


    // Constructors

        //- Construct from components
        inline ReactingParcel
        (
            ReactingCloud<ParcelType>& owner,
            const label typeId,
            const vector& position,
            const label celli,
            const scalar d0,
            const vector& U0,
            const scalar nParticle0,
            const scalarField& YGas0,
            const scalarField& YLiquid0,
            const scalarField& YSolid0,
            const scalarField& YMixture0,
            const constantProperties& constProps
        );

        //- Construct from Istream
        ReactingParcel
        (
            const Cloud<ParcelType>& c,
            Istream& is,
            bool readFields = true
        );

        //- Construct and return a clone
        autoPtr<ParcelType> clone() const
        {
            return autoPtr<ParcelType>(new ReactingParcel<ParcelType>(*this));
        }


    // Member Functions

        // Access

            //- Return mass fractions of gases
            inline const scalarField& YGas() const;
            inline scalarField& YGas();

            //- Return mass fractions of liquids
            inline const scalarField& YLiquid() const;
            inline scalarField& YLiquid();

            //- Return mass fractions of solids
            inline const scalarField& YSolid() const;
            inline scalarField& YSolid();

            //- Return mass fractions of mixture
            inline const scalarField& YMixture() const;
            inline scalarField& YMixture();

            //- Return initial mass
            inline const scalar mass0() const;
            inline scalar& mass0();


        // Tracking

            //- Move parcel
            template<class TrackingData>
            bool move
            (
                TrackingData&
            );


        // I-O

            static void readFields
            (
                ReactingCloud<ParcelType>& c
            );

            static void writeFields
            (
                const ReactingCloud<ParcelType>& c
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ReactingParcelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ReactingParcel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

