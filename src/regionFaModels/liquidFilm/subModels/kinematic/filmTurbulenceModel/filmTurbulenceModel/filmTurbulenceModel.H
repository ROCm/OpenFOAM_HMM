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
    Foam::regionModels::areaSurfaceFilmModels::filmTurbulenceModel

Description
    Base class for film turbulence models

SourceFiles
    filmTurbulenceModel.C
    filmTurbulenceModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef regionFaModels_filmTurbulenceModel_H
#define regionFaModels_filmTurbulenceModel_H

#include "areaFieldsFwd.H"
#include "runTimeSelectionTables.H"
#include "faMatrices.H"
#include "liquidFilmBase.H"
#include "fvMesh.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                       Class filmTurbulenceModel Declaration
\*---------------------------------------------------------------------------*/

class filmTurbulenceModel
{
    // Private Member Functions

        //- No copy construct
        filmTurbulenceModel(const filmTurbulenceModel&) = delete;

        //- No copy assignment
        void operator=(const filmTurbulenceModel&) = delete;


public:

    // Public Enumerations

        //- Options for the friction models
        enum frictionMethodType
        {
            mquadraticProfile,
            mlinearProfile,
            mDarcyWeisbach,
            mManningStrickler
        };

        //- Options for the shear stress models
        enum shearMethodType
        {
            msimple,
            mwallFunction
        };


protected:

    // Protected Data

        //- Reference to liquidFilmBase
        const liquidFilmBase& film_;

        //- Names for friction models
        static const Enum<frictionMethodType> frictionMethodTypeNames_;

        //- Names for shear stress models
        static const Enum<shearMethodType> shearMethodTypeNames_;

        //- Model dictionary
        const dictionary dict_;

        //- Method used
        const frictionMethodType method_;

        //- Shear method used
        const shearMethodType shearMethod_;

        //- Name of density field (optional)
        word rhoName_;

        //- Reference density needed for incompressible calculations
        scalar rhoRef_;


public:

    //- Runtime type information
    TypeName("filmTurbulenceModel");


    // Declare runtime constructor selection table

         declareRunTimeSelectionTable
         (
             autoPtr,
             filmTurbulenceModel,
             dictionary,
             (
                liquidFilmBase& film,
                const dictionary& dict
             ),
             (film, dict)
         );


    // Constructors

        //- Construct from type name, dictionary and surface film model
        filmTurbulenceModel
        (
            const word& modelType,
            liquidFilmBase& film,
            const dictionary& dict
        );


    // Selectors

        //- Return a reference to the selected injection model
        static autoPtr<filmTurbulenceModel> New
        (
            liquidFilmBase& film,
            const dictionary& dict
        );


    //- Destructor
    virtual ~filmTurbulenceModel() = default;


    // Member Functions

        // Access

            //- Return film
            const liquidFilmBase& film() const;


        // Turbulence

            //- Return the effective viscous stress (laminar + turbulent)
            tmp<volSymmTensorField> devRhoReff() const;

            //- Return primary region friction
            tmp<faVectorMatrix> primaryRegionFriction
            (
                areaVectorField& U
            ) const;

            //- Return rho if specified otherwise rhoRef
            tmp<volScalarField> rho() const;

            //- Return the wall film surface friction
            virtual tmp<areaScalarField> Cw() const;

            //- Return the film turbulence viscosity
            virtual tmp<areaScalarField> mut() const = 0;


        // Evaluation

            //- Correct/update the model
            virtual void correct() = 0;

            //- Return the source for the film momentum equation
            virtual tmp<faVectorMatrix> Su(areaVectorField& U) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
