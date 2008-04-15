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
    Foam::ThermoParcel

Description
    Thermodynamic parcel class with one-way coupling with the continuous
    phase.

    Includes Kinematic parcel sub-models, plus:
    - heat transfer

SourceFiles
    thermoParcelI.H
    ThermoParcel.C
    thermoParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef ThermoParcel_H
#define ThermoParcel_H

#include "IOstream.H"
#include "autoPtr.H"
#include "interpolationCellPoint.H"
#include "contiguous.H"

#include "KinematicParcel.H"
#include "ThermoCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ParcelType>
class ThermoParcel;

template<class ParcelType>
Ostream& operator<<
(
    Ostream&,
    const ThermoParcel<ParcelType>&
);

/*---------------------------------------------------------------------------*\
                       Class ThermoParcel Declaration
\*---------------------------------------------------------------------------*/

template <class ParcelType>
class ThermoParcel
:
    public KinematicParcel<ParcelType>
{

    // Private data

        //- Parcel properties

            //- Temperature [K]
            scalar T_;

            //- Specific heat capacity [J/(kg.K)]
            scalar cp_;


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
            const scalar cpc
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
            const scalar cpc
        );

        //- Calculate new particle temperature
        template<class TrackingData>
        scalar calcHeatTransfer
        (
            TrackingData& td,
            const label celli,
            const scalar dt,
            const scalar rhoc,
            const vector& Uc,
            const scalar muc,
            const scalar Tc,
            const scalar cpc,
            scalar& htc
        );


public:

    //- Runtime type information
    TypeName("ThermoParcel");

    friend class Cloud<ParcelType>;


    //- Class to hold particle constant properties
    class constantProperties
    :
        public KinematicParcel<ParcelType>::constantProperties
    {

        // Private data

            //- Particle initial temperature [K]
            const scalar T0_;

            //- Particle specific heat capacity [J/(kg.K)]
            const scalar cp0_;

            //- Particle emissivity [] (radiation)
            const scalar epsilon0_;

            //- Particle scattering factor [] (radiation)
            const scalar f0_;


    public:

        // Constructors
        constantProperties
        (
            const dictionary& dict
        );

        // Member functions

            // Access

            inline const scalar T0() const;
            inline const scalar cp0() const;
            inline const scalar epsilon0() const;
            inline const scalar f0() const;
    };


    //- Class used to pass tracking data to the trackToFace function
    class trackData
    :
        public KinematicParcel<ParcelType>::trackData
    {

        // Private data

            //- Reference to the cloud containing this particle
            ThermoCloud<ParcelType>& cloud_;

            //- Particle constant properties
            const constantProperties& constProps_;

            //- Interpolators for continuous phase fields

                const interpolationCellPoint<scalar>& TInterp_;
                const interpolationCellPoint<scalar>& cpInterp_;


    public:

        // Constructors

            inline trackData
            (
                ThermoCloud<ParcelType>& cloud,
                const constantProperties& constProps,
                const interpolationCellPoint<scalar>& rhoInterp,
                const interpolationCellPoint<vector>& UInterp,
                const interpolationCellPoint<scalar>& muInterp,
                const interpolationCellPoint<scalar>& TInterp,
                const interpolationCellPoint<scalar>& cpInterp,
                const vector& g
            );


        // Member functions

            inline ThermoCloud<ParcelType>& cloud();

            inline const constantProperties& constProps() const;

            inline const interpolationCellPoint<scalar>& TInterp() const;

            inline const interpolationCellPoint<scalar>& cpInterp() const;
    };


    // Constructors

        //- Construct from components
        inline ThermoParcel
        (
            ThermoCloud<ParcelType>& owner,
            const label typeId,
            const vector& position,
            const label celli,
            const scalar d0,
            const vector& U0,
            const scalar nParticle0,
            const constantProperties& constProps
        );

        //- Construct from Istream
        ThermoParcel
        (
            const Cloud<ParcelType>& c,
            Istream& is,
            bool readFields = true
        );

        //- Construct and return a clone
        autoPtr<ParcelType> clone() const
        {
            return autoPtr<ParcelType>(new ThermoParcel<ParcelType>(*this));
        }


    // Member Functions

        // Access

            //- Return temperature
            inline const scalar T() const;
            inline scalar& T();

            //- Return specific heat capacity
            inline const scalar cp() const;
            inline scalar& cp();

            //- Move parcel
            template<class TrackingData>
            bool move
            (
                 TrackingData&
            );


        //- Patch interaction

            //- Overridable function to handle the particle hitting a wallPatch
            template<class TrackingData>
            void hitWallPatch
            (
                 const wallPolyPatch&,
                 TrackingData& td
            );

            void hitWallPatch
            (
                 const wallPolyPatch&,
                 int&
            );


        // I-O

            static void readFields
            (
                ThermoCloud<ParcelType>& c
            );

            static void writeFields
            (
                const ThermoCloud<ParcelType>& c
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ThermoParcelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ThermoParcel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

