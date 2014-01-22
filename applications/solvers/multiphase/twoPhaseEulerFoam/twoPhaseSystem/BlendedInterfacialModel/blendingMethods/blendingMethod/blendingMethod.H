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
    Foam::blendingMethod

Description

SourceFiles
    blendingMethod.C

\*---------------------------------------------------------------------------*/

#ifndef blendingMethod_H
#define blendingMethod_H

#include "dictionary.H"
#include "runTimeSelectionTables.H"
#include "phaseModel.H" 

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class blendingMethod Declaration
\*---------------------------------------------------------------------------*/

class blendingMethod
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        blendingMethod(const blendingMethod&);

        //- Disallow default bitwise assignment
        void operator=(const blendingMethod&);


protected:

    // Protected data

        //- Phase 1
        const phaseModel& phase1_;

        //- Phase 2
        const phaseModel& phase2_;


public:

    //- Runtime type information
    TypeName("blendingMethod");


    // Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        blendingMethod,
        dictionary,
        (
            const dictionary& dict,
            const phaseModel& phase1,
            const phaseModel& phase2
        ),
        (dict, phase1, phase2)
    );


    // Constructors

        //- Construct from a dictionary and two phases
        blendingMethod
        (
            const dictionary& dict,
            const phaseModel& phase1,
            const phaseModel& phase2
        );


    // Selector

        static autoPtr<blendingMethod> New
        (
            const dictionary& dict,
            const phaseModel& phase1,
            const phaseModel& phase2
        );


    //- Destructor
    virtual ~blendingMethod();


    // Member Functions

        //- Factor for first phase
        virtual tmp<volScalarField> f1() const = 0;

        //- Factor for second phase
        virtual tmp<volScalarField> f2() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
