/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 OpenFOAM Foundation
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
    Foam::liftModel

Description

SourceFiles
    liftModel.C
    newLiftModel.C

\*---------------------------------------------------------------------------*/

#ifndef liftModel_H
#define liftModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "dictionary.H"
#include "phaseModel.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class liftModel Declaration
\*---------------------------------------------------------------------------*/

class liftModel
{
protected:

    // Protected data

        const dictionary& dict_;
        const volScalarField& alpha1_;
        const phaseModel& phase1_;
        const phaseModel& phase2_;


public:

    //- Runtime type information
    TypeName("liftModel");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            liftModel,
            dictionary,
            (
                const dictionary& dict,
                const volScalarField& alpha1,
                const phaseModel& phase1,
                const phaseModel& phase2
            ),
            (dict, alpha1, phase1, phase2)
        );


    // Constructors

        //- Construct null
        liftModel();

        //- Construct from components
        liftModel
        (
            const dictionary& dict,
            const volScalarField& alpha1,
            const phaseModel& phase1,
            const phaseModel& phase2
        );


    //- Destructor
    virtual ~liftModel();


    // Selectors

        static autoPtr<liftModel> New
        (
            const dictionary& dict,
            const volScalarField& alpha1,
            const phaseModel& phase1,
            const phaseModel& phase2
        );


    // Member Functions

        //- Lift force
        virtual tmp<volVectorField> F(const volVectorField& U) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
