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
    Foam::injectionModel

Description
    Injection model class for use with surface film modelling

    Sets the newly injected fluid properties by manipulating two arrays:

      - massToInject
      - diameterToInject

    This mass is later subtracted from the film


SourceFiles
    injectionModelI.H
    injectionModel.C

\*---------------------------------------------------------------------------*/

#ifndef injectionModel_H
#define injectionModel_H

#include "surfaceFilmModel.H"
#include "runTimeSelectionTables.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                      Class injectionModel Declaration
\*---------------------------------------------------------------------------*/

class injectionModel
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        injectionModel(const injectionModel&);

        //- Disallow default bitwise assignment
        void operator=(const injectionModel&);


protected:

    // Protected data

        //- Reference to the owner surface film model
        const surfaceFilmModel& owner_;

        //- Model coefficients dictionary
        dictionary coeffs_;


public:

    //- Runtime type information
    TypeName("injectionModel");


    // Declare runtime constructor selection table

         declareRunTimeSelectionTable
         (
             autoPtr,
             injectionModel,
             dictionary,
             (
                const surfaceFilmModel& owner,
                const dictionary& dict
             ),
             (owner, dict)
         );

    // Constructors

        //- Construct null
        injectionModel(const surfaceFilmModel& owner);

        //- Construct from type name, dictionary and surface film model
        injectionModel
        (
            const word& type,
            const surfaceFilmModel& owner,
            const dictionary& dict
        );


    // Selectors

        //- Return a reference to the selected phase change model
        static autoPtr<injectionModel> New
        (
            const surfaceFilmModel& owner,
            const dictionary& dict
        );


    //- Destructor
    virtual ~injectionModel();


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
                scalarField& massToInject,
                scalarField& diameterToInject
            ) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "injectionModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
