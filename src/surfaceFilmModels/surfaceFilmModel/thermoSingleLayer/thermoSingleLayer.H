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
    Foam::thermoSingleLayer

Description
    Thermodynamic form of single-cell layer surface film model

    Note: defining enthalpy as cp(T - Tstd) - when using liquids from the
    thermophysical library, their enthalpies are calculated similarly, where
    Tstd = 298.15 K

SourceFiles
    thermoSingleLayer.C

\*---------------------------------------------------------------------------*/

#ifndef thermoSingleLayer_H
#define thermoSingleLayer_H

#include "kinematicSingleLayer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceFilmModels
{

// Forward declaration of classes
class phaseChangeModel;

/*---------------------------------------------------------------------------*\
                      Class thermoSingleLayer Declaration
\*---------------------------------------------------------------------------*/

class thermoSingleLayer
:
    public kinematicSingleLayer
{
private:

    // Private member functions

        //- Disallow default bitwise copy construct
        thermoSingleLayer(const thermoSingleLayer&);

        //- Disallow default bitwise assignment
        void operator=(const thermoSingleLayer&);


protected:

    // Protected data

        // Thermo properties  TODO: currently read from coeffs dictionary

            //- Specific heat capacity / [J/kg/K]
            volScalarField cp_;


        // Model parameters

            //- Heat transfer coefficient bewteen wall and film [W/m2/K]
            dimensionedScalar htcw_;

            //- Heat transfer coefficient bewteen film surface and primary
            //  region [W/m2/K]
            dimensionedScalar htcs_;


        // Fields

            //- Temperature / [K]
            volScalarField T_;

            //- Sensible)enthalpy / [J/kg]
            volScalarField hs_;


        // Source term fields

            // Film region - registered to the film region mesh
            // Note: need boundary value mapped from primary region, and then
            // pushed into the patch internal field

                //- Energy / [J/m2/s]
                volScalarField hsSp_;


            // Primary region - registered to the primary region mesh
            // Internal use only - not read-in

                //- Energy / [J/m2/s]
                volScalarField hsSpPrimary_;


        // Fields mapped from primary region - registered to the film region
        // Note: need both boundary AND patch internal fields to be mapped

            //- Temperature / [K]
            volScalarField TPrimary_;


        // Sub-models

            //- Phase change
            autoPtr<phaseChangeModel> phaseChange_;


       // Detached surface properties

           //- Energy sink from film to carrier phase [kg]
           scalarField hsSpDetach_;


    // Protected member functions

        //- Initialise the film model - called on construction
        void initialise();

        //- Read control parameters from dictionary
        virtual bool read();

        //- Reset source term fields
        virtual void resetPrimaryRegionSourceTerms();

        //- Transfer fields from the primary region to the film region
        virtual void transferPrimaryRegionFields();

        //- Update the film sub-models
        virtual void updateSubmodels();

        //- Return the wall/surface heat transfer term for the enthalpy equation
        virtual tmp<fvScalarMatrix> q(volScalarField& h) const;


        // Equations

            //- Solve energy equation
            virtual void solveEnergy();


public:

    //- Runtime type information
    TypeName("thermoSingleLayer");


    // Constructors

        //- Construct from components
        thermoSingleLayer
        (
            const word& modelType,
            const fvMesh& mesh,
            const dimensionedVector& g
        );


    //- Destructor
    virtual ~thermoSingleLayer();


    // Member Functions

        // Model parameters

            //- Return the heat transfer coefficient bewteen wall and film
            inline const dimensionedScalar& htcw() const;

            //- Return the Heat transfer coefficient bewteen film surface and
            //  primary region
            inline const dimensionedScalar& htcs() const;


        // Thermo properties

            //- Return sensible enthalpy as a function of temperature
            inline tmp<Foam::volScalarField> hs
            (
                const volScalarField& T
            ) const;

            //- Return temperature as a function of sensible enthalpy
            inline tmp<Foam::volScalarField> T
            (
                const volScalarField& hs
            ) const;


        // Fields

            //- Return the film temperature [K]
            virtual const volScalarField& T() const;

            //- Return the film specific heat capacity [J/kg/K]
            virtual const volScalarField& cp() const;


         // Source fields (read/write access)

            //- External hook to add sources to the film
            virtual void addSources
            (
                const label patchI,            // patchI on primary region
                const label faceI,             // faceI of patchI
                const scalar massSource,       // [kg]
                const vector& momentumSource,  // [kg.m/s] (tangential momentum)
                const scalar pressureSource,   // [kg.m/s] (normal momentum)
                const scalar energySource      // [J]
            );


            // Primary region

                //- Energy / [J/m2/s]
                inline volScalarField& hsSpPrimary();


            // Film region

                //- Energy / [J/m2/s]
                inline volScalarField& hsSp();


        // Sub-models

            //- Phase change
            inline phaseChangeModel& phaseChange();


        // Evolution

            //- Evolve the film equations
            virtual void evolveFilm();


       // I-O

            //- Provide some feedback
            virtual void info() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "thermoSingleLayerI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
