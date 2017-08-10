/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::phaseChangeTwoPhaseMixtures::constantSolidification

Description
    Solidification model. Simple model driven by temperature as:

    mSol = C*rho*alpha*(T - Tsolification)


SourceFiles
    constantSolidification.C

\*---------------------------------------------------------------------------*/

#ifndef constantSolidification_H
#define constantSolidification_H

#include "InterfaceCompositionModel.H"


// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * * * *//

namespace Foam
{
namespace meltingEvaporationModels
{

/*---------------------------------------------------------------------------*\
                           Class constantSolidification
\*---------------------------------------------------------------------------*/

template<class Thermo, class OtherThermo>
class constantSolidification
:
    public InterfaceCompositionModel<Thermo, OtherThermo>
{
    // Private data

        //- Condensation coefficient [1/s/K]
        dimensionedScalar C_;

        //- Phase transition temperature
        const dimensionedScalar Tactivate_;


public:

    //- Runtime type information
    TypeName("constantSolidification");


    // Constructors

        //- Construct from components
        constantSolidification
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~constantSolidification()
    {}


    // Member Functions

        //- Explicit mass transfer coefficient
        virtual tmp<volScalarField> Kexp
        (
            label variable,
            const volScalarField& field
        ) const;

        //- Semi implicit species mass transfer coefficient
        virtual tmp<volScalarField> Kimp
        (
            label variable,
            const volScalarField& field
        ) const;

        //- Return T transition between phases
        virtual const dimensionedScalar& Tactivate() const;

        //- Derivate sign of the source term:
        // d(Sk-i)/d(variable) > 0 or d(Si-k)/d(variable) < 0
        virtual label dSdVariable();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace meltingEvaporationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "constantSolidification.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
