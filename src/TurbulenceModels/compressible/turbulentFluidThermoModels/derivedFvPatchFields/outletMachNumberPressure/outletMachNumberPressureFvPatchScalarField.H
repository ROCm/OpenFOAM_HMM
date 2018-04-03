/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::outletMachNumberPressureFvPatchScalarField

Group
    grpOutletBoundaryConditions

Description
    This BC set the outlet pressure for compressible flows is such as a
    choked condition is achieved given an input Mach number.
    Optionally, if the BC is operated in a non-choked condition (choked false),
    the constants C1 and A1 are used to calculate the corresponding Mach number at
    the outlet.

    The static pressure is calculates as :

        p_s = p_tot / [1+(k-1)/2*M^2]^(k/(k-1))

    where

    p_tot: [Pa] mass averaged total pressure on outlet patch
    k: [-] mass averaged heat capacity ratio on outlet patch
    M: [-] target Mach number on the outlet patch defining as either choked or non-choked

    For choked conditions:

        M = constant defined in dict

    For non-choked conditions:

        M = A1/(C1*A_outlet)*sqrt(2/(k-1)*[r^(2/k)-r^((k+1)/k)], r=p_back/p_tot

    where

        A_outlet: [m2] area of outlet patch
        A1: [m2] constant defined in dict
        C1: [-] constant defined in dict
        p_back: [Pa] constant defined in dict


Usage
    \table
        Property        | Description                 | Required | Default value
        pBack           | Back pressure               | yes
        M               | Average desired mach number | no
        C1              | Model input                 | no       | 0.0
        A1              | Model input                 | no       | 0.0
        choked          | The outlet is considered chocked | yes
        relax           | relaxation factor (1 fully relax) | no | 0.0
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            outletMachNumberPressure;
        pBack           101325;
        c1              1;
        A1              0.008;
        choked          false;
        value           uniform 101325;
    }
    \endverbatim

SourceFiles
    outletMachNumberPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef outletMachNumberPressureFvPatchScalarField_H
#define outletMachNumberPressureFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class plenumPressureFvPatch Declaration
\*---------------------------------------------------------------------------*/

class outletMachNumberPressureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Mach number
        scalar M_;

        //- Back pressure
        scalar pBack_;

        //- Model constant
        scalar c1_;

        //- Model constant area input
        scalar A1_;

        //- The name of the flux field
        word phiName_;

        //- The name of the rho field
        word rhoName_;

        //- The name of the U field
        word UName_;

        //- Choked
        Switch choked_;

        //- Relaxation factor
        scalar relax_;


public:

    //- Runtime type information
    TypeName("outletMachNumberPressure");


    // Constructors

        //- Construct from patch and internal field
        outletMachNumberPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        outletMachNumberPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given outletMachNumberPressureFvPatchScalarField
        //  onto a new patch
        outletMachNumberPressureFvPatchScalarField
        (
            const outletMachNumberPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        outletMachNumberPressureFvPatchScalarField
        (
            const outletMachNumberPressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new outletMachNumberPressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        outletMachNumberPressureFvPatchScalarField
        (
            const outletMachNumberPressureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new outletMachNumberPressureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
