/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2020 OpenCFD Ltd.
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
    Foam::heatTransferCoeffModels::ReynoldsAnalogy

Description
    Heat transfer coefficient calculation based on Reynolds Analogy,
    which is used to relate turbulent momentum and heat transfer.

    The heat transfer coefficient is derived
    from the skin friction coefficient:

    \f[
        C_f = \frac{\tau_w}{0.5 \rho_\infty |U|^2}
    \f]

    as:

    \f[
        h = 0.5 \rho_\infty c_{p,\infty} |U_{\infty}| C_f
    \f]

    where
    \vartable
      h       | Convective heat transfer coefficient of the flow
      \rho_\infty  | Reference fluid density
      c_{p,\infty} | Reference specific heat capacity at constant pressure
      U_{\infty}   | Reference velocity
      C_f          | Skin friction coefficient
      \tau_w       | Wall shear stress
    \endvartable

Usage
    Minimal example by using \c system/controlDict.functions:
    \verbatim
    heatTransferCoeff1
    {
        // Mandatory and other optional entries
        ...
        htcModel        ReynoldsAnalogy;

        // Conditional mandatory entries (runtime modifiable)
        UInf            (10 0 0);

        // Conditional optional entries (runtime modifiable)
        Cp              <CpName>;
        rho             <rhoName>;

        // mandatory if Cp == CpInf
        CpInf           1005;

        // mandatory if rho == rhoInf
        rhoInf          1;
    }
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                            | Type | Reqd | Dflt
      type      | Model name: ReynoldsAnalogy            | word | yes  | -
      UInf      | Reference velocity                   | scalar | yes  | -
      Cp      | Name of reference specific heat capacity | word | no   | Cp
      CpInf   | Reference specific heat capacity value | scalar | cndtnl | -
      rho     | Name of reference fluid density          | word | no   | rho
      rhoInf  | Reference fluid density value          | scalar | cndtnl | -
    \endtable

Note
  - to use a reference \c Cp, set \c Cp to \c CpInf
  - to use a reference \c rho, set \c rho to \c rhoInf

See also
  - Foam::heatTransferCoeffModel
  - Foam::heatTransferCoeffModels::fixedReferenceTemperature
  - Foam::heatTransferCoeffModels::localReferenceTemperature

SourceFiles
    ReynoldsAnalogy.C

\*---------------------------------------------------------------------------*/

#ifndef heatTransferCoeffModels_ReynoldsAnalogy_H
#define heatTransferCoeffModels_ReynoldsAnalogy_H

#include "heatTransferCoeffModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace heatTransferCoeffModels
{

/*---------------------------------------------------------------------------*\
                   Class ReynoldsAnalogy Declaration
\*---------------------------------------------------------------------------*/

class ReynoldsAnalogy
:
    public heatTransferCoeffModel
{
protected:

    // Protected Data

        //- Name of velocity field
        word UName_;

        //- Reference velocity
        vector URef_;

        //- Name of density field
        word rhoName_;

        //- Reference density
        scalar rhoRef_;

        //- Name of specific heat capacity field
        word CpName_;

        //- Reference specific heat capacity
        scalar CpRef_;


    // Protected Member Functions

        virtual tmp<Field<scalar>> rho(const label patchi) const;

        virtual tmp<Field<scalar>> Cp(const label patchi) const;

        virtual tmp<volSymmTensorField> devReff() const;

        tmp<FieldField<Field, scalar>> Cf() const;

        //- Set the heat transfer coefficient
        virtual void htc
        (
            volScalarField& htc,
            const FieldField<Field, scalar>& q
        );


public:

    //- Runtime type information
    TypeName("ReynoldsAnalogy");


    // Constructors

        //- Construct from components
        ReynoldsAnalogy
        (
            const dictionary& dict,
            const fvMesh& mesh,
            const word& TName
        );

        //- No copy construct
        ReynoldsAnalogy(const ReynoldsAnalogy&) = delete;

        //- No copy assignment
        void operator=(const ReynoldsAnalogy&) = delete;


    //- Destructor
    virtual ~ReynoldsAnalogy() = default;


    // Member Functions

        //- Read from dictionary
        virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace heatTransferCoeffModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
