/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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
    Foam::injectionModel

Description
    Base class for film injection models

SourceFiles
    injectionModel.C
    injectionModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef injectionModel_H
#define injectionModel_H

#include "subModelBase.H"
#include "runTimeSelectionTables.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

/*---------------------------------------------------------------------------*\
                      Class injectionModel Declaration
\*---------------------------------------------------------------------------*/

class injectionModel
:
    public subModelBase
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        injectionModel(const injectionModel&);

        //- Disallow default bitwise assignment
        void operator=(const injectionModel&);


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

        //- Return a reference to the selected injection model
        static autoPtr<injectionModel> New
        (
            const surfaceFilmModel& owner,
            const dictionary& dict
        );


    //- Destructor
    virtual ~injectionModel();


    // Member Functions

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
} // End namespace regionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
