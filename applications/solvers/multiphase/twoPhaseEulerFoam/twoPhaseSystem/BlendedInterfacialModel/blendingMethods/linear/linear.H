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
    Foam::linear

Description

SourceFiles
    linear.C

\*---------------------------------------------------------------------------*/

#ifndef linear_H
#define linear_H

#include "blendingMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blendingMethods
{

/*---------------------------------------------------------------------------*\
                         Class linear Declaration
\*---------------------------------------------------------------------------*/

class linear
:
    public blendingMethod
{
    // Private data

        //- Max fraction of phase 1 which can be considered fully dispersed
        const dimensionedScalar maxFullyDispersedAlpha1_;

        //- Max fraction of phase 1 which can be considered partly dispersed
        const dimensionedScalar maxPartlyDispersedAlpha1_;

        //- Max fraction of phase 2 which can be considered fully dispersed
        const dimensionedScalar maxFullyDispersedAlpha2_;

        //- Max fraction of phase 2 which can be considered partly dispersed
        const dimensionedScalar maxPartlyDispersedAlpha2_;


public:

    //- Runtime type information
    TypeName("linear");


    // Constructors

        //- Construct from a dictionary and two phases
        linear
        (
            const dictionary& dict,
            const phaseModel& phase1,
            const phaseModel& phase2
        );


    //- Destructor
    ~linear();


    // Member Functions
        
        //- Factor for primary phase
        virtual tmp<volScalarField> f1() const;

        //- Factor for secondary phase
        virtual tmp<volScalarField> f2() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace blendingMethods
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
