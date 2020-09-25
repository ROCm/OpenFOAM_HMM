/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Foam::uniformFixedValueFaPatchField

Description
    This finiteArea boundary condition
    provides a uniform fixed value condition, which can also be time-varying.

Usage
    \table
        Property     | Description                  | Required | Default
        uniformValue | uniform value                | yes |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            uniformFixedValue;
        uniformValue    constant 0.2;
    }
    \endverbatim

See also
    Foam::Function1Types
    Foam::fixedValueFvPatchField
    Foam::uniformFixedValueFvPatchField.C

SourceFiles
    uniformFixedValueFaPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef uniformFixedValueFaPatchField_H
#define uniformFixedValueFaPatchField_H

#include "fixedValueFaPatchField.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class uniformFixedValueFaPatch Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class uniformFixedValueFaPatchField
:
    public fixedValueFaPatchField<Type>
{
    // Private Data

        //- Value
        autoPtr<Function1<Type>> uniformValue_;


public:

    //- Runtime type information
    TypeName("uniformFixedValue");


    // Constructors

        //- Construct from patch and internal field
        uniformFixedValueFaPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&
        );

        //- Construct from patch, internal field and dictionary
        uniformFixedValueFaPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        uniformFixedValueFaPatchField
        (
            const uniformFixedValueFaPatchField<Type>&,
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const faPatchFieldMapper&
        );

        //- Copy construct
        uniformFixedValueFaPatchField
        (
            const uniformFixedValueFaPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<faPatchField<Type>> clone() const
        {
            return tmp<faPatchField<Type>>
            (
                new uniformFixedValueFaPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        uniformFixedValueFaPatchField
        (
            const uniformFixedValueFaPatchField<Type>&,
            const DimensionedField<Type, areaMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<faPatchField<Type>> clone
        (
            const DimensionedField<Type, areaMesh>& iF
        ) const
        {
            return tmp<faPatchField<Type>>
            (
                new uniformFixedValueFaPatchField<Type>(*this, iF)
            );
        }


    //- Destructor
    virtual ~uniformFixedValueFaPatchField() = default;


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
    #include "uniformFixedValueFaPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
