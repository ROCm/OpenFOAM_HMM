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
    Foam::phaseChangeTwoPhaseMixtures::kineticGasEvaporation

Description
    kineticGasEvaporation evaporation model.

SourceFiles
    kineticGasEvaporation.C

\*---------------------------------------------------------------------------*/

#ifndef kineticGasEvaporation_H
#define kineticGasEvaporation_H


#include "InterfaceCompositionModel.H"
#include "saturationPressureModel.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * * * *//

namespace Foam
{

class phasePair;

namespace meltingEvaporationModels
{

/*---------------------------------------------------------------------------*\
                   Class kineticGasEvaporation
\*---------------------------------------------------------------------------*/

template<class Thermo, class OtherThermo>
class kineticGasEvaporation
:
    public InterfaceCompositionModel<Thermo, OtherThermo>
{
    // Private data

        //- Evaporation coefficient
        dimensionedScalar C_;

        //- Activation temperature
        const dimensionedScalar Tactivate_;

        //- Molar weight of the vapour in the continous phase
        dimensionedScalar Mv_;

        //- Pointer to saturation pressure model
        autoPtr<saturationPressureModel> saturationPressureModel_;

        //- Source
        mutable volScalarField phi0_;

        //- Spread source
        mutable volScalarField sPhi_;

        //- Diffusivity for source term
        dimensionedScalar D_;

    // Private Member Functions

        //- Constant of propotionality between partial pressure and mass
        //  fraction
        tmp<volScalarField> wRatio() const;


public:

    //- Runtime type information
    TypeName("kineticGasEvaporation");


    // Constructors

        //- Construct from components
        kineticGasEvaporation
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~kineticGasEvaporation()
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

        //- Return Tref
        virtual const dimensionedScalar& Tactivate() const;


        //- Explicit mass transfer for energy
        virtual tmp<volScalarField> KexpEnergy
        (
            label modelVariable,
            const volScalarField& field
        ) const;


        //- Derivate sign of the source term:
        // d(Sk-i)/d(variable) > 0 or d(Si-k)/d(variable) < 0
        virtual label dSdVariable();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace meltingEvaporationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "kineticGasEvaporation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif

// ************************************************************************* //
