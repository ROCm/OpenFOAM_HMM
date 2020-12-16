/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::outletMappedUniformInletFvPatchField

Group
    grpInletBoundaryConditions

Description
    The \c outletMappedUniformInlet is an inlet boundary condition that
    - averages the patch field of \<Type\> over a specified "outlet" patch
    and uniformly applies the averaged value over a specified inlet patch.
    - optionally, the averaged value can be scaled
    and/or offset by a pair of specified values.

    The governing equation of the boundary condition is:

    \f[
        \phi_{inlet} = f \phi_{outlet} + \phi_{offset}
    \f]

    where
    \vartable
      \phi_{inlet}   | Spatially-uniform patch-field value at an inlet patch
      \phi_{outlet}  | Averaged patch-field value at an outlet patch
      f              | User-defined fraction value
      \phi_{offset}  | User-defined offset value
    \endvartable

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        // Mandatory entries (unmodifiable)
        type            outletMappedFilterInlet;
        outletPatch     <outletPatchName>;

        // Optional entries (unmodifiable)
        fraction        0.1;
        offset          10;    // (1 0 0);
        phi             phi;

        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property    | Description                          | Type   | Reqd | Dflt
      type        | Type name: outletMappedUniformInlet  | word   | yes  | -
      outletPatch | Name of patch to be mapped           | word   | yes  | -
      fraction    | Fraction value                       | scalar | no   | 1
      offset      | Offset value                         | Type   | no   | Zero
      phi         | Name of operand flux field           | word   | no   | phi
    \endtable

    The inherited entries are elaborated in:
      - \link fixedValueFvPatchFields.H \endlink

See also
  - Foam::fixedValueFvPatchField
  - Foam::outletMappedUniformInletHeatAdditionFvPatchField
  - Foam::outletMappedUniformInletTemperatureFvPatchField

SourceFiles
    outletMappedUniformInletFvPatchField.C
    outletMappedUniformInletFvPatchFields.C

\*---------------------------------------------------------------------------*/

#ifndef outletMappedUniformInletFvPatchField_H
#define outletMappedUniformInletFvPatchField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class outletMappedUniformInletFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class outletMappedUniformInletFvPatchField
:
    public fixedValueFvPatchField<Type>
{
    // Private Data

        //- Name of the outlet patch to be mapped
        word outletPatchName_;

       //- Name of operand flux field
        word phiName_;

        //- Fraction value
        scalar fraction_;

        //- Offset value
        Type offset_;


public:

    //- Runtime type information
    TypeName("outletMappedUniformInlet");


    // Constructors

        //- Construct from patch and internal field
        outletMappedUniformInletFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        outletMappedUniformInletFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given outletMappedUniformInletFvPatchField
        //- onto a new patch
        outletMappedUniformInletFvPatchField
        (
            const outletMappedUniformInletFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        outletMappedUniformInletFvPatchField
        (
            const outletMappedUniformInletFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new outletMappedUniformInletFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        outletMappedUniformInletFvPatchField
        (
            const outletMappedUniformInletFvPatchField<Type>&,
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
                new outletMappedUniformInletFvPatchField<Type>(*this, iF)
            );
        }


    // Member Functions

        // Access

            //- Name of the outlet patch to be mapped
            const word& outletPatchName() const
            {
                return outletPatchName_;
            }


        // Evaluation

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "outletMappedUniformInletFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
