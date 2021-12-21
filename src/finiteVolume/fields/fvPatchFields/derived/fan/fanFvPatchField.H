/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2019 OpenCFD Ltd.
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
    Foam::fanFvPatchField

Group
    grpCoupledBoundaryConditions

Description
    This boundary condition provides a jump condition, using the \c cyclic
    condition as a base.

    The jump is specified as a \c Function1 type, to enable the use of, e.g.
    constant, polynomial, table values.

    The switch nonDimensional can be used for a non-dimensional table,
    in combination with uniformJump = true.
    As inputs it needs the fan RPM (rpm) and the mean diameter (dm).

    The non-dimensional U for the table is calculated as follows:

    \verbatim
        phi = 120*Un/(PI^3*dm*rpm)
        where:
            dm is the mean diameter.
            rpm is the RPM of the fan.
    \endverbatim

    The non-dimensional pressure:

    \verbatim
        Psi = 2 deltaP/(rho*(sqr(PI*omega*dm)))
        where:
            deltaP is the pressure drop
    \endverbatim

    The non-dimensional table should be given as Psi = F(phi).

Usage
    \table
        Property    | Description                           | Required | Default
        patchType   | underlying patch type should be \c cyclic | yes |
        jumpTable   | jump data, e.g. \c csvFile            | yes |
        phi         | flux field name                       | no  | phi
        rho         | density field name                    | no  | rho
        uniformJump | apply uniform pressure based on avg velocity | no | false
        nonDimensional | use non-dimensional table          | no | false
        rpm         | fan rpm (non-dimensional table)       | no |
        dm          | mean diameter (non-dimensional table) | no |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            fan;
        patchType       cyclic;
        jumpTable       csvFile;

        jumpTableCoeffs
        {
            nHeaderLine     1;
            refColumn       0;
            componentColumns 1(1);
            separator       ",";
            mergeSeparators no;
            file            "<constant>/UvsPressure";
        }
        value           uniform 0;
    }
    \endverbatim

    The above example shows the use of a comma separated (CSV) file to specify
    the jump condition.

Note
    The underlying \c patchType should be set to \c cyclic

See also
    Foam::Function1Types

SourceFiles
    fanFvPatchField.C
    fanFvPatchFields.H
    fanFvPatchFields.C
    fanFvPatchFieldsFwd.H

\*---------------------------------------------------------------------------*/

#ifndef fanFvPatchField_H
#define fanFvPatchField_H

#include "uniformJumpFvPatchField.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class fanFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fanFvPatchField
:
    public uniformJumpFvPatchField<Type>
{
    // Private Data

        //- Name of the flux transporting the field
        word phiName_;

        //- Name of the density field for normalising the mass flux if necessary
        word rhoName_;

        //- Apply uniform pressure drop
        bool uniformJump_;

        //- Use non-dimensional curve
        bool nonDimensional_;

        //- Fan rpm (for non-dimensional curve)
        autoPtr<Function1<scalar>> rpm_;

        //- Fan mean diameter (for non-dimensional curve)
        autoPtr<Function1<scalar>> dm_;


    // Private Member Functions

        //- Calculate the fan pressure jump
        void calcFanJump();


public:

    //- Runtime type information
    TypeName("fan");


    // Constructors

        //- Construct from patch and internal field
        fanFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fanFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fanFvPatchField onto a new patch
        fanFvPatchField
        (
            const fanFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fanFvPatchField
        (
            const fanFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new fanFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fanFvPatchField
        (
            const fanFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new fanFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream& os) const;
};


//- Specialisation of the jump-condition for the pressure
template<>
void fanFvPatchField<scalar>::calcFanJump();


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fanFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
