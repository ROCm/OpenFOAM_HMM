/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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
    Foam::kinematicSingleLayer

Description
    Kinematic form of single-cell layer surface film model

SourceFiles
    kinematicSingleLayer.C

\*---------------------------------------------------------------------------*/

#ifndef kinematicSingleLayer_H
#define kinematicSingleLayer_H

#include "surfaceFilmModel.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceFilmModels
{

// Forward declaration of classes
class injectionModel;

/*---------------------------------------------------------------------------*\
                   Class kinematicSingleLayer Declaration
\*---------------------------------------------------------------------------*/

class kinematicSingleLayer
:
    public surfaceFilmModel
{
private:

    // Private member functions

        //- Disallow default bitwise copy construct
        kinematicSingleLayer(const kinematicSingleLayer&);

        //- Disallow default bitwise assignment
        void operator=(const kinematicSingleLayer&);


protected:

    // Protected data

        // Mesh databases

            //- Film region mesh database
            fvMesh filmRegion_;

            //- Patch normal vectors
            volVectorField nHat_;

            //- Face area magnitudes / [m2]
            volScalarField magSf_;

            //- List of patch IDs on the primary region coupled with the film
            //  region
            labelList primaryPatchIDs_;

            //- List of patch IDs on oppositte side of the film region
            labelList filmTopPatchIDs_;

            //- List of patch IDs on the film region coupled with the primary
            //  region
            labelList filmBottomPatchIDs_;


        // Solution parameters

            //- Momentum predictor
            Switch momentumPredictor_;

            //- Number of outer correctors
            label nOuterCorr_;

            //- Number of PISO-like correctors
            label nCorr_;

            //- Number of non-orthogonal correctors
            label nNonOrthCorr_;

            //- Cumulative continuity error
            scalar cumulativeContErr_;


        // Model parameters

            //- Skin frition coefficient for film/primary region interface
            scalar Cf_;

            //- Stable film thickness - film cannot detach until this is reached
            dimensionedScalar deltaStable_;


        // Thermo properties

            // Fields

                //- Initiliased thermo flag
                bool initialisedThermo_;

                //- Density / [kg/m3]
                volScalarField rho_;

                //- Dynamic viscosity / [Pa.s]
                volScalarField mu_;

                //- Surface tension / [m/s2]
                volScalarField sigma_;


        // Fields

            //- Film thickness / [m]
            volScalarField delta_;

            //- Velocity - mean / [m/s]
            volVectorField U_;

            //- Velocity - surface / [m/s]
            volVectorField Us_;

            //- Velocity - wall / [m/s]
            volVectorField Uw_;

            //- Film thickness*density (helper field) / [kg/m2]
            volScalarField deltaRho_;

            //- Mass flux (includes film thickness) / [kg.m/s]
            surfaceScalarField phi_;


            // Transfer fields - to the primary region

                //- Film mass available for transfer
                volScalarField massForPrimary_;

                //- Parcel diameters originating from film
                volScalarField diametersForPrimary_;

                //- Film mass evolved via phase change
                volScalarField massPhaseChangeForPrimary_;


        // Source term fields

            // Film region - registered to the film region mesh
            // Note: need boundary value mapped from primary region, and then
            // pushed into the patch internal field

                //- Momementum / [kg/m/s2]
                volVectorField USp_;

                //- Pressure / [Pa]
                volScalarField pSp_;

                //- Mass / [kg/m2/s]
                volScalarField rhoSp_;


            // Primary region - registered to the primary region mesh
            // Internal use only - not read-in

                //- Momementum / [kg/m/s2]
                volVectorField USpPrimary_;

                //- Pressure / [Pa]
                volScalarField pSpPrimary_;

                //- Mass / [kg/m2/s]
                volScalarField rhoSpPrimary_;


        // Fields mapped from primary region - registered to the film region
        // Note: need both boundary AND patch internal fields to be mapped

            //- Velocity / [m/s]
            volVectorField UPrimary_;

            //- Pressure / [Pa]
            volScalarField pPrimary_;


        // Sub-models

            //- Injection
            autoPtr<injectionModel> injection_;


       // Checks

           //- Cumulative mass added via sources [kg]
           scalar addedMass_;


       // Detached surface properties

           //- Cumulative mass detached [kg]
           scalar detachedMass_;


    // Protected member functions

        //- Initialise the film model - called on construction
        void initialise();

        //- Read control parameters from dictionary
        virtual bool read();

        //- Correct the thermo fields
        virtual void correctThermoFields();

        //- Reset source term fields
        virtual void resetPrimaryRegionSourceTerms();

        //- Transfer fields from the primary region to the film region
        virtual void transferPrimaryRegionFields();

        //- Correct the source terms for film that detaches from film region
        virtual void correctDetachedFilm();

        // Explicit pressure source contribution
        virtual tmp<volScalarField> pu();

        // Implicit pressure source coefficient
        virtual tmp<volScalarField> pp();

        //- Update the film sub-models
        virtual void updateSubmodels();

        //- Courant number evaluation
        virtual scalar CourantNumber() const;

        //- Continuity check
        virtual void continuityCheck();

        //- Update film surface velocities
        virtual void updateSurfaceVelocities();

        //- Return the stress term for the momentum equation
        virtual tmp<fvVectorMatrix> tau(volVectorField& dU) const;

        //- Constrain a film region master/slave boundaries of a field to a
        //  given value
        template<class Type>
        void constrainFilmField
        (
            Type& field,
            const typename Type::cmptType& value
        );


        // Equations

            //- Solve continuity equation
            virtual void solveContinuity();

            //- Solve for film velocity
            virtual tmp<fvVectorMatrix> solveMomentum
            (
                const volScalarField& pu,
                const volScalarField& pp
            );

            //- Solve coupled velocity-thickness equations
            virtual void solveThickness
            (
                const volScalarField& pu,
                const volScalarField& pp,
                const fvVectorMatrix& UEqn
            );


public:

    //- Runtime type information
    TypeName("kinematicSingleLayer");


    // Constructors

        //- Construct from components
        kinematicSingleLayer
        (
            const word& modelType,
            const fvMesh& mesh,
            const dimensionedVector& g
        );


    //- Destructor
    virtual ~kinematicSingleLayer();


    // Member Functions

        // Access

            //- Return the film mesh database
            virtual const fvMesh& film() const;

            //- Return the patch normal vectors
            inline const volVectorField& nHat() const;

            //- Return the face area magnitudes / [m2]
            inline const volScalarField& magSf() const;

            //- Return the list of coupled patches on the film region
            virtual const labelList& filmBottomPatchIDs() const;

            //- Return the list of patches oppositte coupled patches
            virtual const labelList& filmTopPatchIDs() const;

            //- Return the list of coupled patches on the primary region
            virtual const labelList& primaryPatchIDs() const;


        // Solution parameters

            //- Return the momentum predictor
            inline const Switch& momentumPredictor() const;

            //- Return the number of outer correctors
            inline label nOuterCorr() const;

            //- Return the number of PISO correctors
            inline label nCorr() const;

            //- Return the number of non-orthogonal correctors
            inline label nNonOrthCorr() const;


        // Model parameters

            //- Return the skin friction coefficient
            inline scalar Cf() const;


        // Thermo properties

            //- Return const access to the dynamic viscosity / [Pa.s]
            inline const volScalarField& mu() const;

            //- Return const access to the surface tension / [m/s2]
            inline const volScalarField& sigma() const;


        // Fields

            //- Return const access to the film thickness / [m]
            inline const volScalarField& delta() const;

            //- Return the film velocity [m/s]
            virtual const volVectorField& U() const;

            //- Return the film surface velocity [m/s]
            virtual const volVectorField& Us() const;

            //- Return the film wall velocity [m/s]
            virtual const volVectorField& Uw() const;

            //- Return the film density [kg/m3]
            virtual const volScalarField& rho() const;

            //- Return the film mean temperature [K]
            virtual const volScalarField& T() const;

            //- Return the film surface temperature [K]
            virtual const volScalarField& Ts() const;

            //- Return the film wall temperature [K]
            virtual const volScalarField& Tw() const;

            //- Return the film specific heat capacity [J/kg/K]
            virtual const volScalarField& cp() const;

            //- Return the film thermal conductivity [W/m/K]
            virtual const volScalarField& kappa() const;


           // Transfer fields - to the primary region

               //- Return the film mass available for transfer
               virtual const volScalarField& massForPrimary() const;

               //- Return the parcel diameters originating from film
               virtual const volScalarField& diametersForPrimary() const;

               //- Return the film mass evolved via phase change
               virtual const volScalarField& massPhaseChangeForPrimary() const;


        // External helper functions

            //- Return true if patchI on the primary region is a coupled patch
            //  to the film region
            virtual bool isFilmPatch(const label patchI) const;

            //- External hook to add sources to the film
            virtual void addSources
            (
                const label patchI,            // patchI on primary region
                const label faceI,             // faceI of patchI
                const scalar massSource,       // [kg]
                const vector& momentumSource,  // [kg.m/s] (tangential momentum)
                const scalar pressureSource,   // [kg.m/s] (normal momentum)
                const scalar energySource = 0  // [J]
            );


         // Source fields (read/write access)

            // Primary region

                //- Momementum / [kg/m/s2]
                inline volVectorField& USpPrimary();

                //- Pressure / [Pa]
                inline volScalarField& pSpPrimary();

                //- Mass / [kg/m2/s]
                inline volScalarField& rhoSpPrimary();


            // Film region

                //- Momentum / [kg/m/s2]
                inline volVectorField& USp();

                //- Pressure / [Pa]
                inline volScalarField& pSp();

                //- Mass / [kg/m2/s]
                inline volScalarField& rhoSp();

                //- Momentum / [kg/m/s2]
                inline const volVectorField& USp() const;

                //- Pressure / [Pa]
                inline const volScalarField& pSp() const;

                //- Mass / [kg/m2/s]
                inline const volScalarField& rhoSp() const;


        // Fields mapped from primary region

            //- Velocity / [m/s]
            inline const volVectorField& UPrimary() const;

            //- Pressure / [Pa]
            inline const volScalarField& pPrimary() const;


        // Sub-models

            //- Injection
            inline injectionModel& injection();


        // Helper functions

            //- Return the gravity tangential component contributions
            inline tmp<volVectorField> gTan() const;

            //- Return the gravity normal-to-patch component contribution
            inline tmp<volScalarField> gNorm() const;

            //- Return the gravity normal-to-patch component contribution
            //  Clipped so that only non-zero if g & nHat_ < 0
            inline tmp<volScalarField> gNormClipped() const;


        // Evolution

            //- Pre-evolve film hook
            virtual void preEvolveFilm();

            //- Evolve the film equations
            virtual void evolveFilm();


        // Source fields

            // Mapped into primary region

                //- Return total mass source - Eulerian phase only
                virtual tmp<DimensionedField<scalar, volMesh> > Srho() const;

                //- Return mass source for specie i - Eulerian phase only
                virtual tmp<DimensionedField<scalar, volMesh> > Srho
                (
                    const label i
                ) const;

                //- Return enthalpy source - Eulerian phase only
                virtual tmp<DimensionedField<scalar, volMesh> > Sh() const;


       // I-O

            //- Provide some feedback
            virtual void info() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "kinematicSingleLayerTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "kinematicSingleLayerI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
