/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2011 OpenCFD Ltd.
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
    Foam::KinematicCloud

Description
    Templated base class for kinematic cloud

    - holds a 'cloudSolution' class that stores all relevant solution info

    - sub-models:
      - Collision model
      - Dispersion model
      - Drag model
      - Injection model
      - Patch interaction model
      - Post-processing model
      - Surface film model

SourceFiles
    KinematicCloudI.H
    KinematicCloud.C

\*---------------------------------------------------------------------------*/

#ifndef KinematicCloud_H
#define KinematicCloud_H

#include "Cloud.H"
#include "kinematicCloud.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "cachedRandom.H"
#include "fvMesh.H"
#include "volFields.H"
#include "fvMatrices.H"
#include "particleForces.H"

#include "IntegrationSchemesFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class CloudType>
class CollisionModel;

template<class CloudType>
class DispersionModel;

template<class CloudType>
class DragModel;

template<class CloudType>
class InjectionModel;

template<class CloudType>
class PatchInteractionModel;

template<class CloudType>
class PostProcessingModel;

template<class CloudType>
class SurfaceFilmModel;


/*---------------------------------------------------------------------------*\
                       Class KinematicCloud Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class KinematicCloud
:
    public Cloud<ParcelType>,
    public kinematicCloud
{
    // Private data

        //- Cloud copy pointer
        autoPtr<KinematicCloud<ParcelType> > cloudCopyPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        KinematicCloud(const KinematicCloud&);

        //- Disallow default bitwise assignment
        void operator=(const KinematicCloud&);


public:

    // Cloud solution helper class
    class cloudSolution
    {
        // Private Data

            //- Reference to the mesh
            const fvMesh& mesh_;

            //- Dictionary used during construction
            dictionary dict_;

            //- Cloud active flag
            const Switch active_;

            //- Transient flag
            Switch transient_;

            //- Calculation frequency - carrier steps per cloud step
            //  NOTE: Steady operation only
            label calcFrequency_;

            //- Maximum particle Courant number
            //  Max fraction of current cell that can be traversed in a single
            //  step
            scalar maxCo_;

            //- Current cloud iteration
            label iter_;

            //- Cloud solution time step
            scalar deltaT_;


            // Run-time options

                //- Flag to indicate whether parcels are coupled to the carrier
                //  phase, i.e. whether or not to generate source terms for
                //  carrier phase
                Switch coupled_;

                //- Flag to correct cell values with latest transfer information
                //  during the lagrangian timestep
                Switch cellValueSourceCorrection_;

                //- Maximum particle track time [s]
                scalar maxTrackTime_;

                //- Flag to indicate whether coupling source terms should be
                //  reset on start-up/first read
                Switch resetSourcesOnStartup_;


        // Private Member Functions

            //- Disallow default bitwise assignment
            void operator=(const cloudSolution&);


    public:

        // Constructors

            //- Construct null from mesh reference
            cloudSolution(const fvMesh& mesh);

            //- Construct from mesh and dictionary
            cloudSolution(const fvMesh& mesh, const dictionary& dict);

            //- Construct copy
            cloudSolution(const cloudSolution& cs);


        //- Destructor
        virtual ~cloudSolution();


        // Member functions

            //- Read properties from dictionary
            void read();


            // Access

                //- Return relaxation coefficient for field
                scalar relaxCoeff(const word& fieldName) const;

                //- Return semi-implicit flag coefficient for field
                bool semiImplicit(const word& fieldName) const;

                //- Return reference to the mesh
                inline const fvMesh& mesh() const;

                //- Return const access to the dictionary
                inline const dictionary& dict() const;

                //- Return the active flag
                inline const Switch active() const;

                //- Return const access to the transient flag
                inline const Switch transient() const;

                //- Return const access to the steady flag
                inline const Switch steadyState() const;

                //- Return const access to the calculation frequency
                inline label calcFrequency() const;

                //- Return const access to the max particle Courant number
                inline scalar maxCo() const;

                //- Return const access to the current cloud iteration
                inline label iter() const;

                //- Increment and return iter counter
                inline label nextIter();

                //- Return the time step
                inline scalar deltaT() const;

                //- Return const access to the coupled flag
                inline const Switch coupled() const;

                //- Return const access to the cell value correction flag
                inline const Switch cellValueSourceCorrection() const;

                //- Return const access to the particle track time
                inline scalar maxTrackTime() const;

                //- Return const access to the reset sources flag
                inline const Switch resetSourcesOnStartup() const;

                //- Source terms dictionary
                inline const dictionary& sourceTermDict() const;

                //- Interpolation schemes dictionary
                inline const dictionary& interpolationSchemes() const;

                //- Integration schemes dictionary
                inline const dictionary& integrationSchemes() const;


        // Helper functions

            //- Returns true if performing a cloud iteration this calc step
            bool solveThisStep() const;

            //- Returns true if possible to evolve the cloud and sets timestep
            //  parameters
            bool canEvolve();

            //- Returns true if writing this step
            bool output() const;
    };


protected:

    // Protected data

        //- References to the mesh and time databases
        const fvMesh& mesh_;

        //- Dictionary of particle properties
        IOdictionary particleProperties_;

        //- Solution properties
        cloudSolution solution_;

        //- Parcel constant properties
        typename ParcelType::constantProperties constProps_;

        //- Sub-models dictionary
        const dictionary subModelProperties_;

        //- Random number generator - used by some injection routines
        cachedRandom rndGen_;

        //- Cell occupancy information for each parcel, (demand driven)
        autoPtr<List<DynamicList<ParcelType*> > > cellOccupancyPtr_;


        // References to the carrier gas fields

            //- Density [kg/m3]
            const volScalarField& rho_;

            //- Velocity [m/s]
            const volVectorField& U_;

            //- Dynamic viscosity [Pa.s]
            const volScalarField& mu_;


        // Environmental properties

            //- Gravity
            const dimensionedVector& g_;


        //- Optional particle forces
        particleForces forces_;


        // References to the cloud sub-models

            //- Collision model
            autoPtr<CollisionModel<KinematicCloud<ParcelType> > >
                collisionModel_;

            //- Dispersion model
            autoPtr<DispersionModel<KinematicCloud<ParcelType> > >
                dispersionModel_;

            //- Drag transfer model
            autoPtr<DragModel<KinematicCloud<ParcelType> > > dragModel_;

            //- Injector model
            autoPtr<InjectionModel<KinematicCloud<ParcelType> > >
                injectionModel_;

            //- Patch interaction model
            autoPtr<PatchInteractionModel<KinematicCloud<ParcelType> > >
                patchInteractionModel_;

            //- Post-processing model
            autoPtr<PostProcessingModel<KinematicCloud<ParcelType> > >
                postProcessingModel_;

            //- Surface film model
            autoPtr<SurfaceFilmModel<KinematicCloud<ParcelType> > >
                surfaceFilmModel_;


        // Reference to the particle integration schemes

            //- Velocity integration
            autoPtr<vectorIntegrationScheme> UIntegrator_;


        // Sources

            //- Momentum
            autoPtr<DimensionedField<vector, volMesh> > UTrans_;

            //- Coefficient for carrier phase U equation
            autoPtr<DimensionedField<scalar, volMesh> > UCoeff_;


        // Initialisation

            //- Set cloud sub-models
            void setModels();


        // Cloud evolution functions

            //- Solve the cloud - calls all evolution functions
            void solve(typename ParcelType::trackData& td);

            //- Pre-evolve
            void preEvolve();

            //- Build the cellOccupancy
            void buildCellOccupancy();

            //- Update (i.e. build) the cellOccupancy if it has
            //  already been used
            void updateCellOccupancy();

            //- Evolve the cloud
            void evolveCloud(typename ParcelType::trackData& td);

            //- Particle motion
            void motion(typename ParcelType::trackData& td);

            //- Move-collide particles
            void moveCollide(typename ParcelType::trackData& td);

            //- Post-evolve
            void postEvolve();

            //- Reset state of cloud
            void cloudReset(KinematicCloud<ParcelType>& c);


public:

    // Constructors

        //- Construct given carrier gas fields
        KinematicCloud
        (
            const word& cloudName,
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& mu,
            const dimensionedVector& g,
            bool readFields = true
        );

        //- Copy constructor with new name
        KinematicCloud(KinematicCloud<ParcelType>& c, const word& name);

        //- Copy constructor with new name - creates bare cloud
        KinematicCloud
        (
            const fvMesh& mesh,
            const word& name,
            const KinematicCloud<ParcelType>& c
        );

        //- Construct and return clone based on (this) with new name
        virtual autoPtr<Cloud<ParcelType> > clone(const word& name)
        {
            return autoPtr<Cloud<ParcelType> >
            (
                new KinematicCloud(*this, name)
            );
        }

        //- Construct and return bare clone based on (this) with new name
        virtual autoPtr<Cloud<ParcelType> > cloneBare(const word& name) const
        {
            return autoPtr<Cloud<ParcelType> >
            (
                new KinematicCloud(this->mesh(), name, *this)
            );
        }


    //- Destructor
    virtual ~KinematicCloud();


    //- Type of parcel the cloud was instantiated for
    typedef ParcelType parcelType;


    // Member Functions

        // Access

            //- Return a reference to the cloud copy
            inline const KinematicCloud& cloudCopy() const;

            //- If the collision model controls the wall interaction,
            //  then the wall impact distance should be zero.
            //  Otherwise, it should be allowed to be the value from
            //  the Parcel.
            inline bool hasWallImpactDistance() const;


            // References to the mesh and databases

                //- Return refernce to the mesh
                inline const fvMesh& mesh() const;

                //- Return particle properties dictionary
                inline const IOdictionary& particleProperties() const;

                //- Return const access to the solution properties
                inline const cloudSolution& solution() const;

                //- Return access to the solution properties
                inline cloudSolution& solution();

                //- Return the constant properties
                inline const typename ParcelType::constantProperties&
                    constProps() const;

                    //- Return reference to the sub-models dictionary
                inline const dictionary& subModelProperties() const;


            // Cloud data

                //- Return refernce to the random object
                inline cachedRandom& rndGen();

                //- Return the cell occupancy information for each
                //  parcel, non-const access, the caller is
                //  responsible for updating it for its own purposes
                //  if particles are removed or created.
                inline List<DynamicList<ParcelType*> >& cellOccupancy();


            // References to the carrier gas fields

                //- Return carrier gas velocity
                inline const volVectorField& U() const;

                //- Return carrier gas density
                inline const volScalarField& rho() const;

                //- Return carrier gas dynamic viscosity
                inline const volScalarField& mu() const;


            // Environmental properties

                //- Gravity
                inline const dimensionedVector& g() const;


            //- Optional particle forces
            inline const particleForces& forces() const;


            // Sub-models

                //- Return const access to the collision model
                inline const CollisionModel<KinematicCloud<ParcelType> >&
                    collision() const;

                //- Return reference to the collision model
                inline CollisionModel<KinematicCloud<ParcelType> >&
                    collision();

                //- Return const-access to the dispersion model
                inline const DispersionModel<KinematicCloud<ParcelType> >&
                    dispersion() const;

                //- Return reference to the dispersion model
                inline DispersionModel<KinematicCloud<ParcelType> >&
                    dispersion();

                //- Return const-access to the drag model
                inline const DragModel<KinematicCloud<ParcelType> >&
                    drag() const;

                //- Return const access to the injection model
                inline const InjectionModel<KinematicCloud<ParcelType> >&
                    injection() const;

                //- Return reference to the injection model
                inline InjectionModel<KinematicCloud<ParcelType> >&
                    injection();

                //- Return const-access to the patch interaction model
                inline const PatchInteractionModel<KinematicCloud<ParcelType> >&
                    patchInteraction() const;

                //- Return reference to post-processing model
                inline PostProcessingModel<KinematicCloud<ParcelType> >&
                    postProcessing();

                //- Return const-access to the surface film model
                inline const SurfaceFilmModel<KinematicCloud<ParcelType> >&
                    surfaceFilm() const;

                //- Return reference to the surface film model
                inline SurfaceFilmModel<KinematicCloud<ParcelType> >&
                    surfaceFilm();


            // Integration schemes

                //-Return reference to velocity integration
                inline const vectorIntegrationScheme& UIntegrator() const;


            // Sources

                // Momentum

                    //- Return reference to momentum source
                    inline DimensionedField<vector, volMesh>& UTrans();

                    //- Return const reference to momentum source
                    inline const DimensionedField<vector, volMesh>&
                        UTrans() const;

                     //- Return coefficient for carrier phase U equation
                    inline DimensionedField<scalar, volMesh>& UCoeff();

                    //- Return const coefficient for carrier phase U equation
                    inline const DimensionedField<scalar, volMesh>&
                        UCoeff() const;

                    //- Return tmp momentum source term
                    inline tmp<fvVectorMatrix> SU(volVectorField& U) const;


        // Check

            //- Total mass injected
            inline scalar massInjected() const;

            //- Total mass in system
            inline scalar massInSystem() const;

            //- Total linear momentum of the system
            inline vector linearMomentumOfSystem() const;

            //- Total linear kinetic energy in the system
            inline scalar linearKineticEnergyOfSystem() const;

            //- Total rotational kinetic energy in the system
            inline scalar rotationalKineticEnergyOfSystem() const;


            // Fields

                //- Return the particle volume fraction field
                //  Note: for particles belonging to this cloud only
                inline const tmp<volScalarField> theta() const;

                //- Return the particle mass fraction field
                //  Note: for particles belonging to this cloud only
                inline const tmp<volScalarField> alpha() const;

                //- Return the particle effective density field
                //  Note: for particles belonging to this cloud only
                inline const tmp<volScalarField> rhoEff() const;


        // Cloud evolution functions

            //- Check parcel properties
            void checkParcelProperties
            (
                ParcelType& parcel,
                const scalar lagrangianDt,
                const bool fullyDescribed
            );

            //- Store the current cloud state
            void storeState();

            //- Reset the current cloud to the previously stored state
            void restoreState();

            //- Reset the cloud source terms
            void resetSourceTerms();

            //- Relax field
            template<class Type>
            void relax
            (
                DimensionedField<Type, volMesh>& field,
                const DimensionedField<Type, volMesh>& field0,
                const word& name
            ) const;

            //- Apply relaxation to (steady state) cloud sources
            void relaxSources(const KinematicCloud<ParcelType>& cloudOldTime);

            //- Evolve the cloud
            void evolve();

            //- Print cloud information
            void info() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "KinematicCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "KinematicCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
