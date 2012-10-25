/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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
    Foam::uniformTemperatureJumpFvPatchScalarField

Group
    grpCoupledBoundaryConditions

Description
    This boundary condition provides a temperature jump condition across a
    coupled pair of cyclic patches, when solving for energy.

    The jump is specified as a \c DataEntry type, to enable the use of, e.g.
    contant, polynomial, table values.

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        patchType    | underlying patch type should be \c cyclic| yes |
        jumpTable    | jump data, e.g. \c csvFile | yes      |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type            uniformTemperatureJump;
        patchType       cyclic;
        jumpTable       constant 100;
        value           uniform 300;
    }
    \endverbatim

    The above example shows the use of a constant jump condition of 300 K.

Note
     The underlying \c patchType should be set to \c cyclic

SeeAlso
    Foam::uniformJumpFvPatchField
    Foam::energyJumpFvPatchScalarField

SourceFiles
    uniformTemperatureJumpFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef uniformTemperatureJumpFvPatchScalarField_H
#define uniformTemperatureJumpFvPatchScalarField_H

#include "temperatureJumpBase.H"
#include "uniformJumpFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class uniformTemperatureJumpFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class uniformTemperatureJumpFvPatchScalarField
:
    public temperatureJumpBase,
    public uniformJumpFvPatchField<scalar>
{

public:

    //- Runtime type information
    TypeName("uniformTemperatureJump");

    // Constructors

        //- Construct from patch and internal field
        uniformTemperatureJumpFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        uniformTemperatureJumpFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  uniformTemperatureJumpFvPatchScalarField onto a new patch
        uniformTemperatureJumpFvPatchScalarField
        (
            const uniformTemperatureJumpFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        uniformTemperatureJumpFvPatchScalarField
        (
            const uniformTemperatureJumpFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar> > clone() const
        {
            return tmp<fvPatchField<scalar> >
            (
                new uniformTemperatureJumpFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        uniformTemperatureJumpFvPatchScalarField
        (
            const uniformTemperatureJumpFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<scalar> > clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<scalar> >
            (
                new uniformTemperatureJumpFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        virtual tmp<Field<scalar> > jump() const
        {
            return uniformJumpFvPatchField::jump();
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
