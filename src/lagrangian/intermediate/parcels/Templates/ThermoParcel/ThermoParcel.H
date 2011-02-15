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
    Foam::ThermoParcel

Description
    Thermodynamic parcel class with one/two-way coupling with the continuous
    phase. Includes Kinematic parcel sub-models, plus:
    - heat transfer

SourceFiles
    ThermoParcelI.H
    ThermoParcel.C
    ThermoParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef ThermoParcel_H
#define ThermoParcel_H

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

template<class ParcelType>
class ThermoParcel
:
    public KinematicParcel<ParcelType>
{
public:

    //- Class to hold thermo particle constant properties
    class constantProperties
    :
        public KinematicParcel<ParcelType>::constantProperties
    {

        // Private data

            //- Particle initial temperature [K]
            scalar T0_;

            //- Minimum temperature [K]
            scalar TMin_;

            //- Particle specific heat capacity [J/(kg.K)]
            scalar Cp0_;

            //- Particle emissivity [] (radiation)
            scalar epsilon0_;

            //- Particle scattering factor [] (radiation)
            scalar f0_;

            //- Default carrier Prandtl number []
            scalar Pr_;


    public:

        // Constructors

            //- Null constructor
            constantProperties();

            //- Copy constructor
            constantProperties(const constantProperties& cp);

            //- Constructor from dictionary
            constantProperties
            (
                const dictionary& parentDict,
                const bool readFields = true
            );


        // Member functions

            // Access

                //- Return const access to the particle initial temperature [K]
                inline scalar T0() const;

                //- Return const access to minimum temperature [K]
                inline scalar TMin() const;

                //- Return const access to the particle specific heat capacity
                //  [J/(kg.K)]
                inline scalar Cp0() const;

                //- Return const access to the particle emissivity []
                //  Active for radiation only
                inline scalar epsilon0() const;

                //- Return const access to the particle scattering factor []
                //  Active for radiation only
                inline scalar f0() const;

                //- Return const access to the default carrier Prandtl number []
                inline scalar Pr() const;
    };


    //- Class used to pass thermo tracking data to the trackToFace function
    class trackData
    :
        public KinematicParcel<ParcelType>::trackData
    {

        // Private data

            //- Reference to the cloud containing this particle
            ThermoCloud<ParcelType>& cloud_;

            //- Local copy of carrier specific heat field
            //  Cp not stored on acrrier thermo, but returned as tmp<...>
            const volScalarField Cp_;


            // Interpolators for continuous phase fields

                //- Temperature field interpolator
                autoPtr<interpolation<scalar> > TInterp_;

                //- Specific heat capacity field interpolator
                autoPtr<interpolation<scalar> > CpInterp_;

                //- Radiation field interpolator
                autoPtr<interpolation<scalar> > GInterp_;


    public:

        // Constructors

            //- Construct from components
            inline trackData
            (
                ThermoCloud<ParcelType>& cloud,
                typename KinematicParcel<ParcelType>::trackData::trackPart
                    part = KinematicParcel<ParcelType>::trackData::tpLinearTrack
            );


        // Member functions

            //- Return access to the owner cloud
            inline ThermoCloud<ParcelType>& cloud();

            //- Return access to the locally stored carrier Cp field
            inline const volScalarField& Cp() const;

            //- Return const access to the interpolator for continuous
            //  phase temperature field
            inline const interpolation<scalar>& TInterp() const;

            //- Return const access to the interpolator for continuous
            //  phase specific heat capacity field
            inline const interpolation<scalar>& CpInterp() const;

            //- Return const access to the interpolator for continuous
            //  radiation field
            inline const interpolation<scalar>& GInterp() const;
    };


protected:

    // Protected data

        // Parcel properties

            //- Temperature [K]
            scalar T_;

            //- Specific heat capacity [J/(kg.K)]
            scalar Cp_;


        // Cell-based quantities

            //- Temperature [K]
            scalar Tc_;

            //- Specific heat capacity [J/(kg.K)]
            scalar Cpc_;


    // Protected Member Functions

        //- Calculate new particle temperature
        template<class TrackData>
        scalar calcHeatTransfer
        (
            TrackData& td,
            const scalar dt,           // timestep
            const label cellI,         // owner cell
            const scalar Re,           // Reynolds number
            const scalar Pr,           // Prandtl number - surface
            const scalar kappa,        // Thermal conductivity - surface
            const scalar d,            // diameter
            const scalar rho,          // density
            const scalar T,            // temperature
            const scalar Cp,           // specific heat capacity
            const scalar NCpW,         // Sum of N*Cp*W of emission species
            const scalar Sh,           // explicit particle enthalpy source
            scalar& dhsTrans,          // sensible enthalpy transfer to carrier
            scalar& Cuh                // linearised heat transfer coefficient
        );


public:

    // Static data members

        //- String representation of properties
        static string propHeader;

        //- Runtime type information
        TypeName("ThermoParcel");


    friend class Cloud<ParcelType>;


    // Constructors

        //- Construct from owner, position, and cloud owner
        //  Other properties initialised as null
        inline ThermoParcel
        (
            ThermoCloud<ParcelType>& owner,
            const vector& position,
            const label cellI,
            const label tetFaceI,
            const label tetPtI
        );

        //- Construct from components
        inline ThermoParcel
        (
            ThermoCloud<ParcelType>& owner,
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
        ThermoParcel
        (
            const Cloud<ParcelType>& c,
            Istream& is,
            bool readFields = true
        );

        //- Construct as a copy
        ThermoParcel(const ThermoParcel& p);

        //- Construct as a copy
        ThermoParcel(const ThermoParcel& p, const ThermoCloud<ParcelType>& c);

        //- Construct and return a (basic particle) clone
        virtual autoPtr<Particle<ParcelType> > clone() const
        {
            return autoPtr<Particle<ParcelType> >(new ThermoParcel(*this));
        }

        //- Construct and return a (basic particle) clone
        virtual autoPtr<Particle<ParcelType> > clone
        (
            const Cloud<ParcelType>& c
        ) const
        {
            return autoPtr<Particle<ParcelType> >
            (
                new ThermoParcel
                (
                    *this,
                    static_cast<const ThermoCloud<ParcelType>&>(c)
                )
            );
        }


    // Member Functions

        // Access

            //- Return const access to temperature
            inline scalar T() const;

            //- Return const access to specific heat capacity
            inline scalar Cp() const;

            //- Return the parcel sensible enthalpy
            inline scalar hs() const;

            //- Return const access to carrier temperature
            inline scalar Tc() const;

            //- Return const access to carrier specific heat capacity
            inline scalar Cpc() const;


        // Edit

            //- Return access to temperature
            inline scalar& T();

            //- Return access to specific heat capacity
            inline scalar& Cp();


        // Main calculation loop

            //- Set cell values
            template<class TrackData>
            void setCellValues
            (
                TrackData& td,
                const scalar dt,
                const label cellI
            );

            //- Correct cell values using latest transfer information
            template<class TrackData>
            void cellValueSourceCorrection
            (
                TrackData& td,
                const scalar dt,
                const label cellI
            );

            //- Calculate surface thermo properties
            template<class TrackData>
            void calcSurfaceValues
            (
                TrackData& td,
                const label cellI,
                const scalar T,
                scalar& Ts,
                scalar& rhos,
                scalar& mus,
                scalar& Pr,
                scalar& kappas
            ) const;

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
            static void readFields(Cloud<ParcelType>& c);

            //- Write
            static void writeFields(const Cloud<ParcelType>& c);


    // Ostream Operator

        friend Ostream& operator<< <ParcelType>
        (
            Ostream&,
            const ThermoParcel<ParcelType>&
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
