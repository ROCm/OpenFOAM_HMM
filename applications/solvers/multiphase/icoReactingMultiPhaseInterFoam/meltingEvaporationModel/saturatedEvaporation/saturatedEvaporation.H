/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenCFD Ltd.
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
    Foam::phaseChangeTwoPhaseMixtures::saturatedEvaporation

Description
    Evaporation model based on the Saturated species concentration at the
    interface.

    mEvap = C*rho*D*(Ysat - Yi)

SourceFiles
    saturatedEvaporation.C

\*---------------------------------------------------------------------------*/

#ifndef saturatedEvaporation_H
#define saturatedEvaporation_H

#include "saturationPressureModel.H"
#include "InterfaceCompositionModel.H"

// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * * * *//

namespace Foam
{
namespace meltingEvaporationModels
{

/*---------------------------------------------------------------------------*\
                            Class saturatedEvaporation
\*---------------------------------------------------------------------------*/

template<class Thermo, class OtherThermo>
class saturatedEvaporation
:
    public InterfaceCompositionModel<Thermo, OtherThermo>
{
    // Private data

        //- Model coefficient
        dimensionedScalar C_;

        //- Vaporization minimum temperature
        const dimensionedScalar Tactivate_;

        //- Saturated species name in both phases
        word saturatedName_;

        //- Saturated species index in continous phase (receiving mass)
        label continuousSaturatedIndex_;

        //- Saturated species index in  dispersed phase (giving mass)
        label dispersedSaturatedIndex_;

        //- Saturation pressure model
        autoPtr<saturationPressureModel> saturationPressureModel_;


    // Private Member Functions

        //- Constant of propotionality between partial pressure and mass
        //  fraction
        tmp<volScalarField> wRatio() const;


public:

    //- Runtime type information
    TypeName("saturatedEvaporation");


    // Constructors

        //- Construct from components
        saturatedEvaporation(const dictionary& dict, const phasePair& pair);


    //- Destructor
    virtual ~saturatedEvaporation()
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

        //- Interface mass fraction
        virtual tmp<volScalarField> Yf
        (
            const word& speciesName,
            const volScalarField& Tf
        ) const;

        //- Return T transition between phases
        virtual const dimensionedScalar& Tactivate() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace meltingEvaporationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "saturatedEvaporation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
