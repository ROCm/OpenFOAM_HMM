/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2021 OpenCFD Ltd.
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
    Foam::regionModels::liquidFilmBase

Description
    Base class for thermal 2D shells

SourceFiles
    liquidFilmBase.C

\*---------------------------------------------------------------------------*/

#ifndef liquidFilmBase_H
#define liquidFilmBase_H

#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "faCFD.H"
#include "volFieldsFwd.H"
#include "uniformDimensionedFields.H"
#include "regionFaModel.H"
#include "faOptions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                         Class liquidFilmBase Declaration
\*---------------------------------------------------------------------------*/

class liquidFilmBase
:
    public regionFaModel
{
protected:

    // Protected Data

        // Solution parameters

            //- Momentum predictor
            Switch momentumPredictor_;

            //- Number of outer correctors
            label nOuterCorr_;

            //- Number of PISO-like correctors
            label nCorr_;

            //- Number of film thickness correctors
            label nFilmCorr_;

            //- Cumulative continuity error
            scalar cumulativeContErr_;

            //- Smallest numerical thickness
            dimensionedScalar h0_;

            //- Delta wet for sub-models
            dimensionedScalar deltaWet_;


        //- Name of the velocity field
        word UName_;

        //- Name of the pressure field
        word pName_;

        //- Reference absolute pressure
        scalar pRef_;


        // Fields

            //- Film hight
            areaScalarField h_;

            //- Film velocity
            areaVectorField Uf_;

            //- Film pressure
            areaScalarField pf_;

             //- Primary region pressure
            areaScalarField ppf_;

            //- Film momentum flux
            edgeScalarField phif_;

            //- Film height flux
            edgeScalarField phi2s_;

            //- Normal gravity field
            areaScalarField gn_;

            //- Gravity
            uniformDimensionedVectorField g_;


        // Mass exchanage fields from the primary region (lagragian)

            //- Mass
            volScalarField massSource_;

            //- Momentum
            volVectorField momentumSource_;

            //- Normal pressure by particles
            volScalarField pnSource_;

            //- Energy
            volScalarField energySource_;

            //- Total mass added
            scalar addedMassTotal_;


        //- faOptions
        Foam::fa::options& faOptions_;


public:

    //- Runtime type information
    TypeName("liquidFilmBase");


    // Declare runtime constructor selection tables

         declareRunTimeSelectionTable
         (
             autoPtr,
             liquidFilmBase,
             dictionary,
             (
                const word& modelType,
                const fvPatch& patch,
                const dictionary& dict
             ),
             (modelType, patch, dict)
         );


    // Constructors


        //- Construct from type name and mesh and dict
        liquidFilmBase
        (
            const word& modelType,
            const fvPatch& patch,
            const dictionary& dict
        );

        //- No copy construct
        liquidFilmBase(const liquidFilmBase&) = delete;

        //- No copy assignment
        void operator=(const liquidFilmBase&) = delete;


    // Selectors

        //- Return a reference to the selected model using dictionary
        static autoPtr<liquidFilmBase> New
        (
            const fvPatch& patch,
            const dictionary& dict
        );


    //- Destructor
    virtual ~liquidFilmBase();


    // Member Functions

        //- Courant number evaluation
        virtual scalar CourantNumber() const;


        // Helper functions

            //- Wall velocity
            tmp<areaVectorField> Uw() const;

            //- Film surface film velocity
            tmp<areaVectorField> Us() const;

            //- Primary region velocity at film hight. Assume the film to be
            // in the viscous sub-layer
            tmp<areaVectorField> Up() const;

            //- Map primary static pressure
            tmp<areaScalarField> pg() const;

            //- Wet indicator using h0
            tmp<areaScalarField> alpha() const;


        // Access functions

            //- Return faOptions
            Foam::fa::options& faOptions();

            //- Access const reference Uf
            const areaVectorField& Uf() const;

            //- Access const reference gn
            const areaScalarField& gn() const;

            //- Gravity
            const uniformDimensionedVectorField& g() const;

            //- Access const reference h
            const areaScalarField& h() const;

            //- Access to momentum flux
            const edgeScalarField& phif() const;

            //- Access continuity flux
            const edgeScalarField& phi2s() const;


            //- Return h0
            const dimensionedScalar& h0() const;

            //- Access to this region
            const regionFaModel& region() const;

            //- Access to pRef
            scalar pRef();

            //- Name of the U field
            word UName() const;


            // Transfer fields - to the primary region (lagragian injection)

                //- Return mass transfer source - Eulerian phase only
                //virtual tmp<volScalarField> primaryMassTrans() const = 0;

                //- Return the film mass available for transfer to cloud
                virtual const volScalarField& cloudMassTrans() const = 0;

                //- Return the parcel diameters originating from film to cloud
                virtual const volScalarField& cloudDiameterTrans() const = 0;



        // Evolution

            //- Pre-evolve film
            virtual void preEvolveRegion();

            //- Post-evolve film
            virtual void postEvolveRegion();


        // Thermo variables

             //- Access const reference mu
            virtual const areaScalarField& mu() const = 0;

            //- Access const reference rho
            virtual const areaScalarField& rho() const = 0;

            //- Access const reference sigma
            virtual const areaScalarField& sigma() const = 0;

            //- Access const reference Tf
            virtual const areaScalarField& Tf() const = 0;

            //- Access const reference Cp
            virtual const areaScalarField& Cp() const = 0;


        // External hook to add sources and mass exchange variables


            //- Add sources
            virtual void addSources
            (
                const label patchi,            // patchi on primary region
                const label facei,             // facei of patchi
                const scalar massSource,       // [kg]
                const vector& momentumSource,  // [kg.m/s] (tang'l momentum)
                const scalar pressureSource,   // [kg.m/s] (normal momentum)
                const scalar energySource = 0  // [J]
            );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
