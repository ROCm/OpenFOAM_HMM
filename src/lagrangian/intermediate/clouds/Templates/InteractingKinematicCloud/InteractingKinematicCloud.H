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
    Foam::InteractingKinematicCloud

Description
    Templated base class for interactingKinematic cloud

    - InteractingKinematic only
      - Particle collisions
      - Dispersion model
      - Drag model
      - Injection model
      - Wall interaction model

SourceFiles
    InteractingKinematicCloudI.H
    InteractingKinematicCloud.C

\*---------------------------------------------------------------------------*/

#ifndef InteractingKinematicCloud_H
#define InteractingKinematicCloud_H

#include "Cloud.H"
#include "kinematicCloud.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "Random.H"
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
class DispersionModel;

template<class CloudType>
class DragModel;

template<class CloudType>
class InjectionModel;

template<class CloudType>
class CollisionModel;

template<class CloudType>
class PostProcessingModel;

template<class CloudType>
class PatchInteractionModel;

/*---------------------------------------------------------------------------*\
                   Class InteractingKinematicCloud Declaration
\*---------------------------------------------------------------------------*/

template<class ParcelType>
class InteractingKinematicCloud
:
    public Cloud<ParcelType>,
    public kinematicCloud
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        InteractingKinematicCloud(const InteractingKinematicCloud&);

        //- Disallow default bitwise assignment
        void operator=(const InteractingKinematicCloud&);


protected:

    // Protected data

        //- References to the mesh and time databases
        const fvMesh& mesh_;

        //- Dictionary of particle properties
        IOdictionary particleProperties_;

        //- Parcel constant properties
        typename ParcelType::constantProperties constProps_;

        //- Parcel type id - used to flag the type of parcels issued by this
        //  cloud
        const label parcelTypeId_;

        //- Flag to indicate whether parcels are coupled to the carrier phase
        //  i.e. whether or not to generate source terms for carrier phase
        const Switch coupled_;

        //- Flag to correct cell values with latest transfer information
        //  during the lagrangian timestep
        const Switch cellValueSourceCorrection_;

        //- Random number generator - used by some injection routines
        Random rndGen_;


        // References to the carrier gas fields

            //- Density
            const volScalarField& rho_;

            //- Velocity
            const volVectorField& U_;

            //- Dynamic viscosity
            const volScalarField& mu_;


        // Environmental properties

            //- Gravity
            const dimensionedVector& g_;


        //- Optional particle forces
        particleForces forces_;

        //- Interpolation schemes dictionary
        dictionary interpolationSchemes_;


        // References to the cloud sub-models

            //- Dispersion model
            autoPtr<DispersionModel<InteractingKinematicCloud<ParcelType> > >
                dispersionModel_;

            //- Drag transfer model
            autoPtr<DragModel<InteractingKinematicCloud<ParcelType> > >
                dragModel_;

            //- Injector model
            autoPtr<InjectionModel<InteractingKinematicCloud<ParcelType> > >
                injectionModel_;

            //- Collision model
            autoPtr<CollisionModel<InteractingKinematicCloud<ParcelType> > >
                collisionModel_;

            //- Patch interaction model
            autoPtr
            <
                PatchInteractionModel<InteractingKinematicCloud<ParcelType> >
            >
                patchInteractionModel_;

            //- Post-processing model
            autoPtr
            <
                PostProcessingModel<InteractingKinematicCloud<ParcelType> >
            >
                postProcessingModel_;


        // Reference to the particle integration schemes

            //- Velocity integration
            autoPtr<vectorIntegrationScheme> UIntegrator_;


        // Sources

            //- Momentum
            DimensionedField<vector, volMesh> UTrans_;


        // Cloud evolution functions

            //- Pre-evolve
            void preEvolve();

            //- Post-evolve
            void postEvolve();


public:

    // Constructors

        //- Construct given carrier gas fields
        InteractingKinematicCloud
        (
            const word& cloudName,
            const volScalarField& rho,
            const volVectorField& U,
            const volScalarField& mu,
            const dimensionedVector& g,
            bool readFields = true
        );


    //- Destructor
    virtual ~InteractingKinematicCloud();


    //- Type of parcel the cloud was instantiated for
    typedef ParcelType parcelType;


    // Member Functions

        // Access

            // References to the mesh and databases

                //- Return the parcel type id
                inline label parcelTypeId() const;

                //- Return refernce to the mesh
                inline const fvMesh& mesh() const;

                //- Return particle properties dictionary
                inline const IOdictionary& particleProperties() const;

                //- Return the constant properties
                inline const typename ParcelType::constantProperties&
                    constProps() const;


            //- Return coupled flag
            inline const Switch coupled() const;

            //- Return cell value correction flag
            inline const Switch cellValueSourceCorrection() const;

            //- Return refernce to the random object
            inline Random& rndGen();


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


            // Interpolations

                //- Return reference to the interpolation dictionary
                inline const dictionary& interpolationSchemes() const;


            // Sub-models

                //- Return const-access to the dispersion model
                inline
                const DispersionModel<InteractingKinematicCloud<ParcelType> >&
                    dispersion() const;

                //- Return reference to the dispersion model
                inline DispersionModel<InteractingKinematicCloud<ParcelType> >&
                    dispersion();

                //- Return const-access to the drag model
                inline const DragModel<InteractingKinematicCloud<ParcelType> >&
                    drag() const;

                //- Return const access to the injection model
                inline
                const InjectionModel<InteractingKinematicCloud<ParcelType> >&
                    injection() const;

                //- Return reference to the injection model
                inline InjectionModel<InteractingKinematicCloud<ParcelType> >&
                    injection();

                //- Return const access to the collision model
                inline
                const CollisionModel<InteractingKinematicCloud<ParcelType> >&
                    collision() const;

                //- Return reference to the collision model
                inline CollisionModel<InteractingKinematicCloud<ParcelType> >&
                    collision();

                //- Return const-access to the patch interaction model
                inline
                const PatchInteractionModel
                <
                    InteractingKinematicCloud<ParcelType>
                >&
                    patchInteraction() const;

                //- Return reference to post-processing model
                inline
                PostProcessingModel<InteractingKinematicCloud<ParcelType> >&
                    postProcessing();


            // Integration schemes

                //-Return reference to velocity integration
                inline const vectorIntegrationScheme& UIntegrator() const;


            // Sources

                // Momentum

                    //- Return reference to momentum source
                    inline DimensionedField<vector, volMesh>& UTrans();

                    //- Return tmp momentum source term - fully explicit
                    inline tmp<DimensionedField<vector, volMesh> > SU() const;


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

            //- Print cloud information
            void info() const;


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

            //- Reset the spray source terms
            void resetSourceTerms();

            //- Evolve the spray (inject, inject)
            void evolve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "InteractingKinematicCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "InteractingKinematicCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
