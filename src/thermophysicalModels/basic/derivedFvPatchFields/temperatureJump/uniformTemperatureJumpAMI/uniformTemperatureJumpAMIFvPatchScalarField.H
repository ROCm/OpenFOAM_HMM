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
    Foam::uniformTemperatureJumpAMIFvPatchScalarField

Group
    grpCoupledBoundaryConditions

Description
    This boundary condition provides a temperature jump condition across a
    coupled pair of non-conformal cyclic patches using an arbitrary mesh
    interface (AMI), when solving for energy.

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
        type            uniformTemperatureJumpAMI;
        patchType       cyclic;
        jumpTable       constant 100;
        value           uniform 300;
    }
    \endverbatim

    The above example shows the use of a constant jump condition of 300 K.

Note
     The underlying \c patchType should be set to \c cyclicAMI

SeeAlso
    Foam::uniformJumpAMIFvPatchField
    Foam::energyJumpAMIFvPatchScalarField

SourceFiles
    uniformTemperatureJumpAMIFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef uniformTemperatureJumpAMIFvPatchScalarField_H
#define uniformTemperatureJumpAMIFvPatchScalarField_H

#include "temperatureJumpAMIBase.H"
#include "uniformJumpAMIFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class uniformTemperatureJumpAMIFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class uniformTemperatureJumpAMIFvPatchScalarField
:
    public temperatureJumpAMIBase,
    public uniformJumpAMIFvPatchField<scalar>
{

public:

    //- Runtime type information
    TypeName("uniformTemperatureJumpAMI");

    // Constructors

        //- Construct from patch and internal field
        uniformTemperatureJumpAMIFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        uniformTemperatureJumpAMIFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  uniformTemperatureJumpAMIFvPatchScalarField onto a new patch
        uniformTemperatureJumpAMIFvPatchScalarField
        (
            const uniformTemperatureJumpAMIFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        uniformTemperatureJumpAMIFvPatchScalarField
        (
            const uniformTemperatureJumpAMIFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<scalar> > clone() const
        {
            return tmp<fvPatchField<scalar> >
            (
                new uniformTemperatureJumpAMIFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        uniformTemperatureJumpAMIFvPatchScalarField
        (
            const uniformTemperatureJumpAMIFvPatchScalarField&,
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
                new uniformTemperatureJumpAMIFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        virtual tmp<Field<scalar> > jump() const
        {
            return uniformJumpAMIFvPatchField::jump();
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
