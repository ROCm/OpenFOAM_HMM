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
    Foam::phaseChangeModel

Description
    Phase change model for surface film modelling.

SourceFiles
    phaseChangeModelI.H
    phaseChangeModel.C

\*---------------------------------------------------------------------------*/

#ifndef phaseChangeModel_H
#define phaseChangeModel_H

#include "surfaceFilmModel.H"
#include "runTimeSelectionTables.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                      Class phaseChangeModel Declaration
\*---------------------------------------------------------------------------*/

class phaseChangeModel
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        phaseChangeModel(const phaseChangeModel&);

        //- Disallow default bitwise assignment
        void operator=(const phaseChangeModel&);


protected:

    // Protected data

        //- Reference to the owner surface film model
        const surfaceFilmModel& owner_;

        //- Model coefficients dictionary
        dictionary coeffs_;


public:

    //- Runtime type information
    TypeName("phaseChangeModel");


    // Declare runtime constructor selection table

         declareRunTimeSelectionTable
         (
             autoPtr,
             phaseChangeModel,
             dictionary,
             (
                const surfaceFilmModel& owner,
                const dictionary& dict
             ),
             (owner, dict)
         );

    // Constructors

        //- Construct null
        phaseChangeModel(const surfaceFilmModel& owner);

        //- Construct from type name, dictionary and surface film model
        phaseChangeModel
        (
            const word& type,
            const surfaceFilmModel& owner,
            const dictionary& dict
        );


    // Selectors

        //- Return a reference to the selected phase change model
        static autoPtr<phaseChangeModel> New
        (
            const surfaceFilmModel& owner,
            const dictionary& dict
        );


    //- Destructor
    virtual ~phaseChangeModel();


    // Member Functions

        // Access

            //- Return the reference to the owner surface film model
            inline const surfaceFilmModel& owner() const;

            //- Return the model coefficients dictionary
            inline const dictionary& coeffs() const;


        // Evolution

            //- Correct
            virtual void correct
            (
                const scalar dt,
                scalarField& dMass,
                scalarField& dEnergy
            ) = 0;


        // Input/output

            //- Output model statistics
            virtual void info() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "phaseChangeModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
