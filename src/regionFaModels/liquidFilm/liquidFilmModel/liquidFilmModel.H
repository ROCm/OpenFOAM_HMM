/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::regionFaModels::liquidFilmModel

Description
    Thin Model Film model.


SourceFiles
    liquidFilmModel.C
    kinematicThinFilmI.H

\*---------------------------------------------------------------------------*/

#ifndef liquidFilmModel_H
#define liquidFilmModel_H

#include "volFieldsFwd.H"
#include "liquidFilmBase.H"
#include "faMesh.H"
#include "filmTurbulenceModel.H"
#include "liquidMixtureProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                        Class liquidFilmModel Declaration
\*---------------------------------------------------------------------------*/

class liquidFilmModel
:
    public liquidFilmBase
{
private:

    // Private member functions

        //- No copy construct
        liquidFilmModel(const liquidFilmModel&) = delete;

        //- No copy assignment
        void operator=(const liquidFilmModel&) = delete;

        //- Initialize liquidFilmModel
        void init();


protected:


    // Thermo properties


            //- Liquid thermo
            liquidMixtureProperties thermo_;


            // Fields

                //- Density [kg/m3]
                areaScalarField rho_;

                //- Dynamic viscosity [Pa.s]
                areaScalarField mu_;

                //- Surface tension [m/s2]
                areaScalarField sigma_;

                //- Film rho*height
                areaScalarField hRho_;

                // Mass exchange sources

                    //- Mass source
                    areaScalarField rhoSp_;

                    //- Momentum source
                    areaVectorField USp_;

                    //- Normal pressure by particles
                    areaScalarField pnSp_;


    // General properties



            //- Turbulence model
            autoPtr<filmTurbulenceModel> turbulence_;


    // Protected member functions


        //- Read control parameters from dictionary
        virtual bool read(const dictionary& dict);


public:

    //- Runtime type information
    TypeName("liquidFilmModel");


    // Constructors


        //- Construct from components and dict
        liquidFilmModel
        (
            const word& modelType,
            const fvPatch& patch,
            const dictionary& dict
        );


    //- Destructor
    virtual ~liquidFilmModel();


    // Member Functions

        // Help function


            //- Correct thermo
            void correctThermoFields();


        // Access functions

            //- Access const reference mu
            const areaScalarField& mu() const;

            //- Access const reference rho
            const areaScalarField& rho() const;

            //- Access const reference sigma
            const areaScalarField& sigma() const;


        // Evolution

            //- Pre-evolve film
            virtual void preEvolveRegion();

            //- Evolve the film
            //virtual void evolveRegion();

            //- Post-evolve film
            virtual void postEvolveRegion();


       // I-O

            //- Provide some feedback
            //virtual void info();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
