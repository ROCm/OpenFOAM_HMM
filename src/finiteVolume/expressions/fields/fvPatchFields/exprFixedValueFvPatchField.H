/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
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
    Foam::exprFixedValueFvPatchField

Description
    A fixed value boundary condition with expressions.

Usage
    \table
        Property     | Description                          | Required | Default
        value        | fixed value                          | yes |
        valueExpr    | expression for uniformValue          | yes |
    \endtable

Note
    Can also just use uniformFixedValueFvPatchField with an expression
    for the PatchFunction1.

SourceFiles
    exprFixedValueFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef exprFixedValueFvPatchField_H
#define exprFixedValueFvPatchField_H

#include "fixedValueFvPatchField.H"
#include "patchExprFieldBase.H"
#include "patchExprDriver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class exprFixedValueFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class exprFixedValueFvPatchField
:
    public fixedValueFvPatchField<Type>,
    public expressions::patchExprFieldBase
{
    //- The parent boundary condition type
    typedef fixedValueFvPatchField<Type> parent_bctype;


protected:

    // Protected Data

        //- Dictionary contents for the boundary condition
        dictionary dict_;

        //- The expression driver
        expressions::patchExpr::parseDriver driver_;


    // Protected Member Functions

        //- Set debug ON if "debug" is enabled
        void setDebug();


public:

    //- Runtime type information
    TypeName("exprFixedValue");


    // Constructors

        //- Construct from patch and internal field
        exprFixedValueFvPatchField
        (
            const fvPatch& p,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        exprFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary& dict,
            const bool valueRequired=true
        );

        //- Construct by mapping onto a new patch
        exprFixedValueFvPatchField
        (
            const exprFixedValueFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        exprFixedValueFvPatchField
        (
            const exprFixedValueFvPatchField<Type>&
        );


        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new exprFixedValueFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        exprFixedValueFvPatchField
        (
            const exprFixedValueFvPatchField<Type>&,
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
                new exprFixedValueFvPatchField<Type>(*this, iF)
            );
        }


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "exprFixedValueFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
