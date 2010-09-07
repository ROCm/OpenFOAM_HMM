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
    Foam::surfaceFilmModel

Description
    Top level class that provides surface film modelling functionality

SourceFiles
    surfaceFilmModelI.H
    surfaceFilmModel.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceFilmModel_H
#define surfaceFilmModel_H

#include "IOdictionary.H"
#include "Switch.H"
#include "dimensionedVector.H"
#include "runTimeSelectionTables.H"
#include "volFieldsFwd.H"
#include "DimensionedField.H"
#include "labelList.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;
class Time;

namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                      Class surfaceFilmModel Declaration
\*---------------------------------------------------------------------------*/

class surfaceFilmModel
:
    public IOdictionary
{
public:

    // Data types

        //- Enumeration listing the possible thermo types
        enum thermoModelType
        {
            tmConstant,
            tmSingleComponent
        };

        //- Named enumeration for the thermoType
        static const NamedEnum<thermoModelType, 2> thermoModelTypeNames_;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        surfaceFilmModel(const surfaceFilmModel&);

        //- Disallow default bitwise assignment
        void operator=(const surfaceFilmModel&);


protected:

    // Protected data

        //- Reference to the mesh database
        const fvMesh& mesh_;

        //- Reference to the time database
        const Time& time_;

        //- Active flag
        Switch active_;

        //- Acceleration due to gravity [m/s2]
        const dimensionedVector& g_;

        //- Name of film region
        word filmRegionName_;

        //- Model coefficients dictionary
        dictionary coeffs_;

        //- Thermo type
        thermoModelType thermoModel_;


    // Protected member functions

        //- Read control parameters from dictionary
        virtual bool read();


public:

    //- Runtime type information
    TypeName("surfaceFilmModel");


    // Declare runtime constructor selection table

         declareRunTimeSelectionTable
         (
             autoPtr,
             surfaceFilmModel,
             mesh,
             (
                const word& type,
                const fvMesh& mesh,
                const dimensionedVector& g
             ),
             (type, mesh, g)
         );

    // Constructors

        //- Construct null
        surfaceFilmModel(const fvMesh& mesh, const dimensionedVector& g);

        //- Construct from type name and fvMesh
        surfaceFilmModel
        (
            const word& type,
            const fvMesh& mesh,
            const dimensionedVector& g
        );


    // Selectors

        //- Return a reference to the selected surface film model
        static autoPtr<surfaceFilmModel> New
        (
            const fvMesh& mesh,
            const dimensionedVector& g
        );


    //- Destructor
    virtual ~surfaceFilmModel();


    // Member Functions

        // Access

            //- Return the reference to the primary mesh database
            const fvMesh& mesh() const;

            //- Return the reference to the time database
            const Time& time() const;

            //- Return the active flag
            inline const Switch& active() const;

            //- Return the model coefficients dictionary
            inline const dictionary& coeffs() const;

            //- Return the thermo type
            inline const thermoModelType& thermoModel() const;

            //- Return the film mesh database
            virtual const fvMesh& film() const = 0;

            //- Return the list of coupled patches on the film region
            virtual const labelList& filmBottomPatchIDs() const = 0;

            //- Return the list of patches oppositte coupled patches
            virtual const labelList& filmTopPatchIDs() const = 0;

            //- Return the list of coupled patches on the primary region
            virtual const labelList& primaryPatchIDs() const = 0;

            //- Return true if primary patchI is coupled to the film region
            virtual bool isFilmPatch(const label patchI) const = 0;

            //- External hook to add sources to the film
            virtual void addSources
            (
                const label patchI,
                const label faceI,
                const scalar massSource,
                const vector& momentumSource,
                const scalar pressureSource,
                const scalar energySource
            ) = 0;


       // Fields

           //- Return the film velocity [m/s]
           virtual const volVectorField& U() const = 0;

           //- Return the film surface velocity [m/s]
           virtual const volVectorField& Us() const = 0;

           //- Return the film wall velocity [m/s]
           virtual const volVectorField& Uw() const = 0;

           //- Return the film density [kg/m3]
           virtual const volScalarField& rho() const = 0;

           //- Return the film mean temperature [K]
           virtual const volScalarField& T() const = 0;

           //- Return the film surface temperature [K]
           virtual const volScalarField& Ts() const = 0;

           //- Return the film wall temperature [K]
           virtual const volScalarField& Tw() const = 0;

            //- Return the film specific heat capacity [J/kg/K]
           virtual const volScalarField& cp() const = 0;

           //- Return the film thermal conductivity [W/m/K]
           virtual const volScalarField& kappa() const = 0;


           // Transfer fields - to the primary region

               //- Return the film mass available for transfer
               virtual const volScalarField& massForPrimary() const = 0;

               //- Return the parcel diameters originating from film
               virtual const volScalarField& diametersForPrimary() const = 0;

               //- Return the film mass evolved via phase change
               virtual const volScalarField& massPhaseChangeForPrimary()
                   const = 0;


       // Evolution

            //- Evolve the film
            virtual void evolve();

            //- Pre-evolve film hook
            virtual void preEvolveFilm() = 0;

            //- Evolve the film
            virtual void evolveFilm() = 0;


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
            virtual void info() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "surfaceFilmModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
