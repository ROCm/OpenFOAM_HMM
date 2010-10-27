/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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
    Foam::KinematicParcel

Description
    Kinematic parcel class with rotational motion (as spherical
    particles only) and one/two-way coupling with the continuous
    phase.

    Sub-models include:
    - drag
    - turbulent dispersion
    - wall interactions
    - many-body collisions, including memory of data from previous collision

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

#include "KinematicCloud.H"
#include "CollisionRecordList.H"
#include "labelFieldIOField.H"
#include "vectorFieldIOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

typedef CollisionRecordList<vector, vector> collisionRecordList;
typedef vectorFieldCompactIOField pairDataFieldCompactIOField;
typedef vectorFieldCompactIOField wallDataFieldCompactIOField;

template<class ParcelType>
class KinematicParcel;

// Forward declaration of friend functions

template<class ParcelType>
Ostream& operator<<
(
    Ostream&,
    const KinematicParcel<ParcelType>&
);

/*---------------------------------------------------------------------------*\
                         Class KinematicParcel Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class KinematicParcel
:
    public Particle<ParcelType>
{
public:

    //- Class to hold kinematic particle constant properties
    class constantProperties
    {
        // Private data

            //- Constant properties dictionary
            const dictionary dict_;

            //- Parcel type id - used for post-processing to flag the type
            //  of parcels issued by this cloud
            const label parcelTypeId_;

            //- Minimum density [kg/m3]
            const scalar rhoMin_;

            //- Particle density [kg/m3] (constant)
            const scalar rho0_;

            //- Minimum particle mass [kg]
            const scalar minParticleMass_;

            //- Young's modulus [N/m2]
            const scalar youngsModulus_;

            //- Poisson's ratio
            const scalar poissonsRatio_;


    public:

        //- Constructor
        constantProperties(const dictionary& parentDict);

        // Member functions

            //- Return const access to the constant properties dictionary
            inline const dictionary& dict() const;

            //- Return const access to the parcel type id
            inline label parcelTypeId() const;

            //- Return const access to the minimum density
            inline scalar rhoMin() const;

            //- Return const access to the particle density
            inline scalar rho0() const;

            //- Return const access to the minimum particle mass
            inline scalar minParticleMass() const;

            //- Return const access to Young's Modulus
            inline scalar youngsModulus() const;

            //- Return const access to Poisson's ratio
            inline scalar poissonsRatio() const;
    };


    //- Class used to pass kinematic tracking data to the trackToFace function
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
            KinematicCloud<ParcelType>& cloud_;


            // Interpolators for continuous phase fields

                //- Density interpolator
                autoPtr<interpolation<scalar> > rhoInterp_;

                //- Velocity interpolator
                autoPtr<interpolation<vector> > UInterp_;

                //- Dynamic viscosity interpolator
                autoPtr<interpolation<scalar> > muInterp_;


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
                KinematicCloud<ParcelType>& cloud,
                trackPart part = tpLinearTrack
            );


        // Member functions

            //- Return access to the owner cloud
            inline KinematicCloud<ParcelType>& cloud();

            //- Return conat access to the interpolator for continuous
            //  phase density field
            inline const interpolation<scalar>& rhoInterp() const;

            //- Return conat access to the interpolator for continuous
            //  phase velocity field
            inline const interpolation<vector>& UInterp() const;

            //- Return conat access to the interpolator for continuous
            //  phase dynamic viscosity field
            inline const interpolation<scalar>& muInterp() const;

            // Return const access to the gravitational acceleration vector
            inline const vector& g() const;

            //- Return the part of the tracking operation taking place
            inline trackPart part() const;

            //- Return access to the part of the tracking operation taking place
            inline trackPart& part();
    };


protected:

    // Protected data

        // Parcel properties

            //- Active flag - tracking inactive when active = false
            bool active_;

            //- Parcel type id
            label typeId_;

            //- Number of particles in Parcel
            scalar nParticle_;

            //- Diameter [m]
            scalar d_;

            //- Target diameter [m]
            scalar dTarget_;

            //- Velocity of Parcel [m/s]
            vector U_;

            //- Force on particle due to collisions [N]
            vector f_;

            //- Angular momentum of Parcel in global reference frame
            // [kg m2/s]
            vector angularMomentum_;

            //- Torque on particle due to collisions in global
            //  reference frame [Nm]
            vector torque_;

            //- Density [kg/m3]
            scalar rho_;

            //- Age [s]
            scalar age_;

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


    // Protected Member Functions

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
            vector& dUTrans,           // momentum transfer to carrier
            scalar& Cud                // linearised drag coefficient
        ) const;


public:

    // Static data members

        //- String representation of properties
        static string propHeader;

        //- Runtime type information
        TypeName("KinematicParcel");


    friend class Cloud<ParcelType>;


    // Constructors

        //- Construct from owner, position, and cloud owner
        //  Other properties initialised as null
        inline KinematicParcel
        (
            KinematicCloud<ParcelType>& owner,
            const vector& position,
            const label cellI,
            const label tetFaceI,
            const label tetPtI
        );

        //- Construct from components
        inline KinematicParcel
        (
            KinematicCloud<ParcelType>& owner,
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
        KinematicParcel
        (
            const Cloud<ParcelType>& c,
            Istream& is,
            bool readFields = true
        );

        //- Construct as a copy
        KinematicParcel(const KinematicParcel& p);

        //- Construct as a copy
        KinematicParcel
        (
            const KinematicParcel& p,
            const KinematicCloud<ParcelType>& c
        );

        //- Construct and return a (basic particle) clone
        virtual autoPtr<Particle<ParcelType> > clone() const
        {
            return autoPtr<Particle<ParcelType> >(new KinematicParcel(*this));
        }

        //- Construct and return a (basic particle) clone
        virtual autoPtr<Particle<ParcelType> > clone
        (
            const Cloud<ParcelType>& c
        ) const
        {
            return autoPtr<Particle<ParcelType> >
            (
                new KinematicParcel
                (
                    *this,
                    static_cast<const KinematicCloud<ParcelType>&>(c)
                )
            );
        }


    // Member Functions

        // Access

            //- Return const access to active flag
            inline bool active() const;

            //- Return const access to type id
            inline label typeId() const;

            //- Return const access to number of particles
            inline scalar nParticle() const;

            //- Return const access to diameter
            inline scalar d() const;

            //- Return const access to target diameter
            inline scalar dTarget() const;

            //- Return const access to velocity
            inline const vector& U() const;

            //- Return const access to force
            inline const vector& f() const;

            //- Return const access to angular momentum
            inline const vector& angularMomentum() const;

            //- Return const access to torque
            inline const vector& torque() const;

            //- Return const access to density
            inline scalar rho() const;

            //- Return const access to the age
            inline scalar age() const;

            //- Return const access to time spent in turbulent eddy
            inline scalar tTurb() const;

            //- Return const access to turbulent velocity fluctuation
            inline const vector& UTurb() const;

            //- Return const access to the collision records
            inline const collisionRecordList& collisionRecords() const;


        // Edit

            //- Return const access to active flag
            inline bool& active();

            //- Return access to type id
            inline label& typeId();

            //- Return access to number of particles
            inline scalar& nParticle();

            //- Return access to diameter
            inline scalar& d();

            //- Return access to target diameter
            inline scalar& dTarget();

            //- Return access to velocity
            inline vector& U();

            //- Return access to force
            inline vector& f();

            //- Return access to angular momentum
            inline vector& angularMomentum();

            //- Return access to torque
            inline vector& torque();

            //- Return access to density
            inline scalar& rho();

            //- Return access to the age
            inline scalar& age();

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

            //- Particle volume
            inline scalar volume() const;

            //- Particle volume for a given diameter
            inline static scalar volume(const scalar d);

            //- Particle projected area
            inline scalar areaP() const;

            //- Projected area for given diameter
            inline static scalar areaP(const scalar d);

            //- Particle surface area
            inline scalar areaS() const;

            //- Surface area for given diameter
            inline static scalar areaS(const scalar d);

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
            bool move(TrackData& td, const scalar trackTime);


        // Patch interactions

            //- Overridable function to handle the particle hitting a face
            //  without trackData
            void hitFace(int& td);

            //- Overridable function to handle the particle hitting a face
            template<class TrackData>
            void hitFace(TrackData& td);

            //- Overridable function to handle the particle hitting a patch
            //  Executed before other patch-hitting functions
            template<class TrackData>
            bool hitPatch
            (
                const polyPatch& p,
                TrackData& td,
                const label patchI,
                const scalar trackFraction,
                const tetIndices& tetIs
            );

            //- Overridable function to handle the particle hitting a patch
            //  Executed before other patch-hitting functions without trackData
            bool hitPatch
            (
                const polyPatch& p,
                int& td,
                const label patchI,
                const scalar trackFraction,
                const tetIndices& tetIs
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
                TrackData& td,
                const tetIndices&
            );

            //- Overridable function to handle the particle hitting a wallPatch
            //  without trackData
            void hitWallPatch
            (
                const wallPolyPatch&,
                int&,
                const tetIndices&
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
