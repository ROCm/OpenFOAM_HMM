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
    Foam::KinematicParcel

Description
    Kinematic parcel class with one-way coupling with the continuous
    phase.

    Sub-models include:
    - drag
    - break-up
    - wall interactions

SourceFiles
    KinematicParcelI.H
    KinematicParcel.C
    KinematicParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef KinematicParcel_H
#define KinematicParcel_H

#include "Particle.H"
#include "IOstream.H"
#include "autoPtr.H"
#include "interpolationCellPoint.H"
#include "contiguous.H"
#include "mathematicalConstants.H"

#include "KinematicCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ParcelType>
class KinematicParcel;

template<class ParcelType>
Ostream& operator<<
(
    Ostream&,
    const KinematicParcel<ParcelType>&
);

/*---------------------------------------------------------------------------*\
                         Class KinematicParcel Declaration
\*---------------------------------------------------------------------------*/

template <class ParcelType>
class KinematicParcel
:
    public Particle<ParcelType>
{

    // Private member data

        //- Parcel type id
        label typeId_;

        //- Diameter [m]
        scalar d_;

        //- Velocity of Parcel [m/s]
        vector U_;

        //- Relative velocity of Parcel [m/s]
        vector Ur_;

        //- Number of particles in Parcel
        scalar nParticle_;

        //- Density [kg/m3]
        scalar rho_;

        //- Time spent in turbulent eddy
        scalar tTurb_;

        //- Turbulent velocity fluctuation
        vector UTurb_;


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
            const scalar muc
        );

        template<class TrackingData>
        void calcUncoupled
        (
            TrackingData& td,
            const scalar dt,
            const scalar rhoc,
            vector& Uc,
            const scalar muc
        );

        //- Calculate new particle velocity
        template<class TrackingData>
        vector calcVelocity
        (
            TrackingData& td,
            const scalar dt,
            const scalar rhoc,
            vector& Uc,
            const scalar muc,
            scalar& Cud
        );


public:

    //- Runtime type information
    TypeName("KinematicParcel");

    friend class Cloud<ParcelType>;


    //- Class to hold particle constant properties
    class constantProperties
    {

        // Private data

            //- Particle density [kg/m3] (constant)
            const scalar rho0_;

            //- Minimum particle mass [kg]
            const scalar minParticleMass_;


    public:

        // Constructors
        constantProperties
        (
            const dictionary& dict
        );

        // Member functions

            // Access

            inline const scalar rho0() const;

            inline const scalar minParticleMass() const;
    };


    //- Class used to pass tracking data to the trackToFace function
    class trackData
    :
        public Particle<ParcelType>::trackData
    {

        // Private data

            //- Reference to the cloud containing this particle
            KinematicCloud<ParcelType>& cloud_;

            //- Particle constant properties
            const constantProperties& constProps_;

            //- Interpolators for continuous phase fields

                const interpolationCellPoint<scalar>& rhoInterp_;
                const interpolationCellPoint<vector>& UInterp_;
                const interpolationCellPoint<scalar>& muInterp_;

            //- Local gravitational or other body-force acceleration
            const vector& g_;


   public:

        bool switchProcessor;
        bool keepParticle;


        // Constructors

            inline trackData
            (
                KinematicCloud<ParcelType>& cloud,
                const constantProperties& constProps,
                const interpolationCellPoint<scalar>& rhoInterp,
                const interpolationCellPoint<vector>& UInterp,
                const interpolationCellPoint<scalar>& muInterp,
                const vector& g
            );


        // Member functions

            inline KinematicCloud<ParcelType>& cloud();

            inline const constantProperties& constProps() const;

            inline const interpolationCellPoint<scalar>& rhoInterp() const;

            inline const interpolationCellPoint<vector>& UInterp() const;

            inline const interpolationCellPoint<scalar>& muInterp() const;

            inline const vector& g() const;
    };


    // Constructors

        //- Construct from components
        inline KinematicParcel
        (
            KinematicCloud<ParcelType>& owner,
            const label typeId,
            const vector& position,
            const label celli,
            const scalar d0,
            const vector& U0,
            const scalar nParticle0,
            const constantProperties& constProps
        );

        //- Construct from Istream
        KinematicParcel
        (
            const Cloud<ParcelType>& c,
            Istream& is,
            bool readFields = true
        );

        //- Construct and return a clone
        autoPtr<ParcelType> clone() const
        {
            return autoPtr<ParcelType>(new KinematicParcel<ParcelType>(*this));
        }


    // Member Functions

        // Access

            //- Return type id
            inline const label& typeId() const;

            //- Return diameter
            inline const scalar& d() const;
            inline scalar& d();

            //- Return velocity
            inline const vector& U() const;
            inline vector& U();

            //- Return relative velocity
            inline const vector& Ur() const;
            inline vector& Ur();

            //- Return number of particles
            inline const scalar& nParticle() const;
            inline scalar& nParticle();

            //- Return density
            inline const scalar& rho() const;
            inline scalar& rho();

            //- Return time spent in turbulent eddy
            inline const scalar& tTurb() const;
            inline scalar& tTurb();

            //- Return turbulent velocity fluctuation
            inline const vector& UTurb() const;
            inline vector& UTurb();

            //- The nearest distance to a wall that
            //  the particle can be in the n direction
            inline scalar wallImpactDistance(const vector& n) const;

            //- Return the index of the face to be used in the interpolation routine
            inline label faceInterpolation() const;

            //- Particle volume
            inline scalar volume() const;

            //- Particle mass
            inline scalar mass() const;

            //- Particle projected area
            inline scalar areaP() const;

            //- Particle surface area
            inline scalar areaS() const;


        // Tracking

            //- Move the parcel
            template<class TrackingData>
            bool move
            (
                TrackingData& td
            );

            //- Overridable function to handle the particle hitting a
            //  processorPatch
            template<class TrackingData>
            void hitProcessorPatch
            (
                const processorPolyPatch&,
                TrackingData& td
            );

            //- Overridable function to handle the particle hitting a
            //  processorPatch without trackData
            virtual void hitProcessorPatch
            (
                const processorPolyPatch&,
                int&
            );

            //- Overridable function to handle the particle hitting a wallPatch
            template<class TrackingData>
            void hitWallPatch
            (
                const wallPolyPatch&,
                TrackingData& td
            );

            //- Overridable function to handle the particle hitting a wallPatch
            //  without trackData
            virtual void hitWallPatch
            (
                const wallPolyPatch&,
                int&
            );

            //- Overridable function to handle the particle hitting a polyPatch
            template<class TrackingData>
            void hitPatch
            (
                const polyPatch&,
                TrackingData& td
            );

            //- Overridable function to handle the particle hitting a polyPatch
            //- without trackData
            void hitPatch
            (
                const polyPatch&,
                int&
            );

            //- Transform the physical properties of the particle
            //  according to the given transformation tensor
            void transformProperties(const tensor& T);

            //- Transform the physical properties of the particle
            //  according to the given separation vector
            void transformProperties(const vector& separation);


        // I-O

            static void readFields
            (
                KinematicCloud<ParcelType>& c
            );

            static void writeFields
            (
                const KinematicCloud<ParcelType>& c
            );

    // Ostream Operator

        friend Ostream& operator<< <ParcelType>
        (
            Ostream&,
            const KinematicParcel<ParcelType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "KinematicParcelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineParcelTypeNameAndDebug(Type, DebugSwitch)                       \
    template<>                                                                \
    const Foam::word KinematicParcel<Type>::typeName(#Type);                  \
    template<>                                                                \
    int KinematicParcel<Type>::debug                                          \
    (                                                                         \
        Foam::debug::debugSwitch(#Type, DebugSwitch)                          \
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "KinematicParcel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

