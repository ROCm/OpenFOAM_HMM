/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::meltingEvaporationModels::kineticGasEvaporation

Description
    Considering the Hertz Knudsen formula, which gives the
    evaporation-condensation flux based on the kinetic theory for flat
    interface:

    \f[
        Flux = C \sqrt{\frac{M}{2 \pi R T_{activate}}}(p - p_{sat})
    \f]

        where:
        \vartable
            Flux    | mass flux rate [kg/s/m2]
            M       | molecular weight
            T_{activate} | saturation temperature
            C       | accomodation coefficient
            R       | universal gas constant
            p_{sat} | saturation pressure
            p       | vapor partial pressure
        \endvartable

    The Clapeyron-Clausius equation relates the pressure to the temperature
    for the saturation condition:

    \f[
        \frac{dp}{dT} = - \frac{L}{T (\nu_v - \nu_l)}
    \f]

        where:
        \vartable
            L       | latent heat
            \nu_v   | inverse of the vapor density
            \nu_l   | inverse of the liquid density
         \endvartable


    Using the above relations:

    \f[
        Flux =
           2 \frac{C}{2 - C}
           \sqrt{\frac{M}{2 \pi R T_{activate}}}
           L (\rho_{v}*\rho_{l}/(\rho_{l} - \rho_{v}))
           (T - T_{activate})/T_{activate}
    \f]

    This assumes liquid and vapour are in equilibrium, then the accomodation
    coefficient are equivalent for the interface. This relation is known as the
    Hertz-Knudsen-Schrage.

    Based on the reference:
    - Van P. Carey, Liquid-Vapor Phase Change Phenomena, ISBN 0-89116836,
      1992, pp. 112-121.


Usage

    Example usage:
    \verbatim
    massTransferModel
    (
        (liquid to gas)
        {
            type                kineticGasEvaporation;
            species             vapour.gas;
            C                   0.1;
            alphaMin            0.0;
            alphaMax            0.2;
            Tactivate           373;
        }
    );
    \endverbatim

    where:
    \table
        Property    | Description             | Required    | Default value
        C           | Accomodation coefficient (C > 0 for evaporation, C < 0 for
            condensation) | yes
        alphaMin    | Minimum value of alpha  | no          | 0.0
        alphaMax    | Maximum values of alpha | no          | 0.5
        Tactivate   | Saturation temperature  | yes
        species     | Specie name on the other phase | yes
    \endtable

SourceFiles
    kineticGasEvaporation.C

\*---------------------------------------------------------------------------*/

#ifndef meltingEvaporationModels_kineticGasEvaporation_H
#define meltingEvaporationModels_kineticGasEvaporation_H


#include "InterfaceCompositionModel.H"

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

        //- 'To' phase maximum value for the mass transfer
        scalar alphaMax_;

        //- 'To' phase minumum value for the mass transfer
        scalar alphaMin_;

        //- Alpha maximum for the rest of phases
        scalar alphaRestMax_;


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
    virtual ~kineticGasEvaporation() = default;


    // Member Functions

        //- Explicit mass transfer coefficient
        virtual tmp<volScalarField> Kexp
        (
            label variable,
            const volScalarField& field
        );

        //- Return Tactivate
        virtual const dimensionedScalar& Tactivate() const;
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
