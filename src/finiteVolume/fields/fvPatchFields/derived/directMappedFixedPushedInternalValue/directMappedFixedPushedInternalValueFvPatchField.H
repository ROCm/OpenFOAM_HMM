/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::directMappedFixedPushedInternalValueFvPatchField

Description
    Recycles the boundary values of a neighbour patch field to the boundary
    and internal values of *this.

SourceFiles
    directMappedFixedPushedInternalValueFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef directMappedFixedPushedInternalValueFvPatchField_H
#define directMappedFixedPushedInternalValueFvPatchField_H

#include "directMappedFixedValueFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class directMappedFixedPushedInternalValueFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class directMappedFixedPushedInternalValueFvPatchField
:
    public directMappedFixedValueFvPatchField<Type>
{

public:

    //- Runtime type information
    TypeName("directMappedFixedPushedInternalValue");


    // Constructors

        //- Construct from patch and internal field
        directMappedFixedPushedInternalValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        directMappedFixedPushedInternalValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given a
        //  directMappedFixedPushedInternalValueFvPatchField  onto a new patch
        directMappedFixedPushedInternalValueFvPatchField
        (
            const directMappedFixedPushedInternalValueFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        directMappedFixedPushedInternalValueFvPatchField
        (
            const directMappedFixedPushedInternalValueFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new directMappedFixedPushedInternalValueFvPatchField<Type>
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        directMappedFixedPushedInternalValueFvPatchField
        (
            const directMappedFixedPushedInternalValueFvPatchField<Type>&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type> > clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type> >
            (
                new directMappedFixedPushedInternalValueFvPatchField<Type>
                (
                    *this,
                    iF
                )
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

#ifdef NoRepository
#   include "directMappedFixedPushedInternalValueFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
