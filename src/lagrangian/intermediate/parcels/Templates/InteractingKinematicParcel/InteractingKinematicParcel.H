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
    Foam::InteractingKinematicParcel

Description
    InteractingKinematic parcel class with one/two-way coupling with
    the continuous phase.

    Sub-models include:
    - drag
    - turbulent dispersion
    - wall interactions

    Includes many-body collisions and rotational motion.  Only
    spherical particles currently supported.

SourceFiles
    InteractingKinematicParcelI.H
    InteractingKinematicParcel.C
    InteractingKinematicParcelIO.C

\*---------------------------------------------------------------------------*/

#ifndef InteractingKinematicParcel_H
#define InteractingKinematicParcel_H

#include "Particle.H"
#include "IOstream.H"
#include "autoPtr.H"
#include "interpolationCellPoint.H"
#include "contiguous.H"

#include "InteractingKinematicCloud.H"
#include "collisionRecordList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class ParcelType>
class InteractingKinematicParcel;

// Forward declaration of friend functions

template<class ParcelType>
Ostream& operator<<
(
    Ostream&,
    const InteractingKinematicParcel<ParcelType>&
);

/*---------------------------------------------------------------------------*\
                  Class InteractingKinematicParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class InteractingKinematicParcel
:
    public Particle<ParcelType>
{
public:

    //- Class to hold interactingKinematic particle constant properties
    class constantProperties
    {
        // Private data

            //- Constant properties dictionary
            const dictionary dict_;

            //- Minimum density [kg/m3]
            const scalar rhoMin_;

            //- Particle density [kg/m3] (constant)
            const scalar rho0_;

            //- Minimum particle mass [kg]
            const scalar minParticleMass_;


    public:

        //- Constructor
        constantProperties(const dictionary& parentDict);

        // Member functions

            //- Return const access to the constant properties dictionary
            inline const dictionary& dict() const;

            //- Return const access to the minimum density
            inline scalar rhoMin() const;

            //- Return const access to the particle density
            inline scalar rho0() const;

            //- Return const access to the minimum particle mass
            inline scalar minParticleMass() const;
    };


    //- Class used to pass interactingKinematic tracking data to the
    //  trackToFace function
    class trackData
    :
        public Particle<ParcelType>::trackData
    {
    public:

        enum trackPart
        {
            tpVelocityHalfStep,
            tpLinearTrack,
            tpRotationalTrack
        };


    private:

        // Private data

            //- Reference to the cloud containing this particle
            InteractingKinematicCloud<ParcelType>& cloud_;

            //- Particle constant properties
            const constantProperties& constProps_;

            // Interpolators for continuous phase fields

                //- Density interpolator
                const interpolation<scalar>& rhoInterp_;

                //- Velocity interpolator
                const interpolation<vector>& UInterp_;

                //- Dynamic viscosity interpolator
                const interpolation<scalar>& muInterp_;

            //- Local gravitational or other body-force acceleration
            const vector& g_;

            // label specifying which part of the integration
            // algorithm is taking place
            trackPart part_;


   public:

        // Constructors

           //- Construct from components
           inline trackData
            (
                InteractingKinematicCloud<ParcelType>& cloud,
                const constantProperties& constProps,
                const interpolation<scalar>& rhoInterp,
                const interpolation<vector>& UInterp,
                const interpolation<scalar>& muInterp,
                const vector& g,
                trackPart part = tpLinearTrack
            );


        // Member functions

            //- Return access to the owner cloud
            inline InteractingKinematicCloud<ParcelType>& cloud();

            //- Return const access to the constant properties
            inline const constantProperties& constProps() const;

            //- Return conat access to the interpolator for continuous
            //  phase density field
            inline const interpolation<scalar>& rhoInterp() const;

            //- Return conat access to the interpolator for continuous
            //  phase velocity field
            inline const interpolation<vector>& UInterp() const;

            //- Return conat access to the interpolator for continuous
            //  phase dynamic viscosity field
            inline const interpolation<scalar>& muInterp() const;

            //- Return const access to the gravitational acceleration vector
            inline const vector& g() const;

            //- Return the part of the tracking operation taking place
            inline trackPart part() const;

            //- Return access to the part of the tracking operation taking place
            inline trackPart& part();
    };


protected:

    // Protected data

        // Parcel properties

            //- Parcel type id
            label typeId_;

            //- Number of particles in Parcel
            scalar nParticle_;

            //- Diameter [m]
            scalar d_;

            //- Velocity of Parcel [m/s]
            vector U_;

            //- Force on particle due to collisions [N]
            vector f_;

            //- Angular momentum of Parcel in global reference frame
            // [kg m^2 s-1]
            vector pi_;

            //- Torque on particle due to collisions in global
            //  reference frame [Nm]
            vector tau_;

            //- Density [kg/m3]
            scalar rho_;

            //- Time spent in turbulent eddy [s]
            scalar tTurb_;

            //- Turbulent velocity fluctuation [m/s]
            vector UTurb_;

            //- Particle collision records
            collisionRecordList collisionRecords_;


        // Cell-based quantities

            //- Density [kg/m3]
            scalar rhoc_;

            //- Velocity [m/s]
            vector Uc_;

            //- Viscosity [Pa.s]
            scalar muc_;


    // Protected member functions

        //- Calculate new particle velocity
        template<class TrackData>
        const vector calcVelocity
        (
            TrackData& td,
            const scalar dt,           // timestep
            const label cellI,         // owner cell
            const scalar Re,           // Reynolds number
            const scalar mu,           // local carrier viscosity
            const scalar d,            // diameter
            const vector& U,           // velocity
            const scalar rho,          // density
            const scalar mass,         // mass
            const vector& Su,          // explicit particle momentum source
            vector& dUTrans            // momentum transfer to carrier
        ) const;


public:

    // Static data members

        //- String representation of properties
        static string propHeader;

        //- Runtime type information
        TypeName("InteractingKinematicParcel");


    friend class Cloud<ParcelType>;


    // Constructors

        //- Construct from owner, position, and cloud owner
        //  Other properties initialised as null
        inline InteractingKinematicParcel
        (
            InteractingKinematicCloud<ParcelType>& owner,
            const vector& position,
            const label cellI
        );

        //- Construct from components
        inline InteractingKinematicParcel
        (
            InteractingKinematicCloud<ParcelType>& owner,
            const vector& position,
            const label cellI,
            const label typeId,
            const scalar nParticle0,
            const scalar d0,
            const vector& U0,
            const vector& f0,
            const vector& pi0,
            const vector& tau0,
            const constantProperties& constProps
        );

        //- Construct from Istream
        InteractingKinematicParcel
        (
            const Cloud<ParcelType>& c,
            Istream& is,
            bool readFields = true
        );

        //- Construct as a copy
        InteractingKinematicParcel(const InteractingKinematicParcel& p);

        //- Construct and return a clone
        autoPtr<InteractingKinematicParcel> clone() const
        {
            return autoPtr<InteractingKinematicParcel>
            (
                new InteractingKinematicParcel(*this)
            );
        }


    // Member Functions

        // Access

            //- Return const access to type id
            inline label typeId() const;

            //- Return const access to number of particles
            inline scalar nParticle() const;

            //- Return const access to diameter
            inline scalar d() const;

            //- Return const access to velocity
            inline const vector& U() const;

            //- Return const access to force
            inline const vector& f() const;

            //- Return const access to angular momentum
            inline const vector& pi() const;

            //- Return const access to torque
            inline const vector& tau() const;

            //- Return const access to density
            inline scalar rho() const;

            //- Return const access to time spent in turbulent eddy
            inline scalar tTurb() const;

            //- Return const access to turbulent velocity fluctuation
            inline const vector& UTurb() const;

            //- Return const access to the collision records
            inline const collisionRecordList& collisionRecords() const;


        // Edit

            //- Return access to type id
            inline label typeId();

            //- Return access to number of particles
            inline scalar& nParticle();

            //- Return access to diameter
            inline scalar& d();

            //- Return access to velocity
            inline vector& U();

            //- Return access to force
            inline vector& f();

            //- Return access to angular momentum
            inline vector& pi();

            //- Return access to torque
            inline vector& tau();

            //- Return access to density
            inline scalar& rho();

            //- Return access to time spent in turbulent eddy
            inline scalar& tTurb();

            //- Return access to turbulent velocity fluctuation
            inline vector& UTurb();

            //- Return access to collision records
            inline collisionRecordList& collisionRecords();


        // Helper functions

            //- The nearest distance to a wall that
            //  the particle can be in the n direction
            inline scalar wallImpactDistance(const vector& n) const;

            //- Return the index of the face to be used in the interpolation
            //  routine
            inline label faceInterpolation() const;

            //- Cell owner mass
            inline scalar massCell(const label cellI) const;

            //- Particle mass
            inline scalar mass() const;

            //- Particle moment of inertia around diameter axis
            inline scalar momentOfInertia() const;

            //- Particle angular velocity
            inline vector omega() const;

            //- Particle radius
            inline scalar r() const;

            //- Particle volume
            inline scalar volume() const;

            //- Particle volume for a given diameter
            inline scalar volume(const scalar d) const;

            //- Particle projected area
            inline scalar areaP() const;

            //- Projected area for given diameter
            inline scalar areaP(const scalar d) const;

            //- Particle surface area
            inline scalar areaS() const;

            //- Surface area for given diameter
            inline scalar areaS(const scalar d) const;

            //- Reynolds number
            inline scalar Re
            (
                const vector& U,        // particle velocity
                const scalar d,         // particle diameter
                const scalar rhoc,      // carrier density
                const scalar muc        // carrier dynamic viscosity
            ) const;


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

            //- Update parcel properties over the time interval
            template<class TrackData>
            void calc
            (
                TrackData& td,
                const scalar dt,
                const label cellI
            );


        // Tracking

            //- Move the parcel
            template<class TrackData>
            bool move(TrackData& td);


        // Patch interactions

            //- Overridable function to handle the particle hitting a patch
            //  Executed before other patch-hitting functions
            template<class TrackData>
            bool hitPatch
            (
                const polyPatch& p,
                TrackData& td,
                const label patchI
            );


            //- Overridable function to handle the particle hitting a patch
            //  Executed before other patch-hitting functions without trackData
            bool hitPatch
            (
                const polyPatch& p,
                int& td,
                const label patchI
            );


            //- Overridable function to handle the particle hitting a
            //  processorPatch
            template<class TrackData>
            void hitProcessorPatch
            (
                const processorPolyPatch&,
                TrackData& td
            );

            //- Overridable function to handle the particle hitting a
            //  processorPatch without trackData
            void hitProcessorPatch
            (
                const processorPolyPatch&,
                int&
            );

            //- Overridable function to handle the particle hitting a wallPatch
            template<class TrackData>
            void hitWallPatch
            (
                const wallPolyPatch&,
                TrackData& td
            );

            //- Overridable function to handle the particle hitting a wallPatch
            //  without trackData
            void hitWallPatch
            (
                const wallPolyPatch&,
                int&
            );

            //- Overridable function to handle the particle hitting a polyPatch
            template<class TrackData>
            void hitPatch
            (
                const polyPatch&,
                TrackData& td
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

            //- Read
            static void readFields(Cloud<ParcelType>& c);

            //- Write
            static void writeFields(const Cloud<ParcelType>& c);


    // Ostream Operator

        friend Ostream& operator<< <ParcelType>
        (
            Ostream&,
            const InteractingKinematicParcel<ParcelType>&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "InteractingKinematicParcelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define defineParcelTypeNameAndDebug(Type, DebugSwitch)                       \
    template<>                                                                \
    const Foam::word InteractingKinematicParcel<Type>::typeName(#Type);       \
    template<>                                                                \
    int InteractingKinematicParcel<Type>::debug                               \
    (                                                                         \
        Foam::debug::debugSwitch(#Type, DebugSwitch)                          \
    );

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "InteractingKinematicParcel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
