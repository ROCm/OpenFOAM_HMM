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
    Foam::meltingEvaporationModels::Lee

Description
    Mass tranfer Lee model. Simple model driven by field value difference as:

    \f[
        \dot{m} = C \rho \alpha (T - T_{activate})/T_{activate}
    \f]

    where C is a model constant.

    if C > 0:
    \f[
        \dot{m} = C \rho \alpha (T - T_{activate})/T_{activate}
    \f]
        for \f[ T > T_{activate} \f]

    and

        \f[ mDot = 0.0 \f]  for \f[  T < T_{activate} \f]


    if C < 0:
    \f[
        \dot{m} = -C \rho \alpha (T_{activate} - T)/T_{activate}
    \f]
        for \f[ T < T_{activate} \f]

    and
        \f[ \dot{m} = 0.0 \f] for \f[ T > T_{activate} \f]

    Based on the reference:
    -# W. H. Lee. "A Pressure Iteration Scheme for Two-Phase Modeling".
    Technical Report LA-UR 79-975. Los Alamos Scientific Laboratory,
    Los Alamos, New Mexico. 1979.

Usage
    Example usage:
    \verbatim
        massTransferModel
        (
            (solid to liquid)
            {
                type            Lee;
                C               40;
                Tactivate       302.78;
            }
        );
    \endverbatim

    Where:

    \table
        Property    | Description            | Required    | Default value
        Tactivate   | Activation temperature | yes
        C           | Model constant         | yes
    \endtable

SourceFiles
    Lee.C

\*---------------------------------------------------------------------------*/

#ifndef meltingEvaporationModels_Lee_H
#define meltingEvaporationModels_Lee_H

#include "InterfaceCompositionModel.H"


// * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * * * * * *//

namespace Foam
{
namespace meltingEvaporationModels
{

/*---------------------------------------------------------------------------*\
                                  Class Lee
\*---------------------------------------------------------------------------*/

template<class Thermo, class OtherThermo>
class Lee
:
    public InterfaceCompositionModel<Thermo, OtherThermo>
{
    // Private data

        //- Condensation coefficient [1/s]
        dimensionedScalar C_;

        //- Phase transition temperature
        const dimensionedScalar Tactivate_;

        //- Phase minimum value for activation
        scalar alphaMin_;


public:

    //- Runtime type information
    TypeName("Lee");


    // Constructors

        //- Construct from components
        Lee
        (
            const dictionary& dict,
            const phasePair& pair
        );


    //- Destructor
    virtual ~Lee() = default;


    // Member Functions

        //- Explicit mass transfer coefficient
        virtual tmp<volScalarField> Kexp
        (
            label variable,
            const volScalarField& field
        );

        //- Return T transition between phases
        virtual const dimensionedScalar& Tactivate() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace meltingEvaporationModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "Lee.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
