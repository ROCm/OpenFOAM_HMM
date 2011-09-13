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
    Foam::selfContainedMappedFixedValueFvPatchField

Description
    Self-contained version of mapped. Does not use information on
    patch, instead holds it locally (and possibly duplicate) so use
    normal mapped in preference and only use this if you cannot
    change the underlying patch type to mapped.

SourceFiles
    selfContainedMappedFixedValueFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef selfContainedMappedFixedValueFvPatchField_H
#define selfContainedMappedFixedValueFvPatchField_H

#include "mappedPatchBase.H"
#include "fixedValueFvPatchFields.H"
#include "interpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class selfContainedMappedFixedValueFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class selfContainedMappedFixedValueFvPatchField
:
    public mappedPatchBase,
    public fixedValueFvPatchField<Type>
{
    // Private data

        //- Name of field to sample - defaults to field associated with this
        //  patchField if not specified
        word fieldName_;

        //- If true adjust the mapped field to maintain average value average_
        const bool setAverage_;

        //- Average value the mapped field is adjusted to maintain if
        //  setAverage_ is set true
        const Type average_;

        //- Interpolation scheme to use for nearestcell mode
        word interpolationScheme_;

        mutable autoPtr<interpolation<Type> > interpolator_;


    // Private Member Functions

        //- Field to sample. Either on my or nbr mesh
        const GeometricField<Type, fvPatchField, volMesh>& sampleField() const;

        //- Access the interpolation method
        const interpolation<Type>& interpolator() const;


public:

    //- Runtime type information
    TypeName("selfContainedMapped");


    // Constructors

        //- Construct from patch and internal field
        selfContainedMappedFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        selfContainedMappedFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct from patch, internal field and distance for normal type
        //  sampling
        selfContainedMappedFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,

            // mappedPatchBase
            const word& sampleRegion,
            const sampleMode sampleMode,
            const word& samplePatch,
            const scalar distance,

            // My settings
            const word& fieldName,
            const bool setAverage,
            const Type average,
            const word& interpolationScheme
        );

        //- Construct by mapping given
        //  selfContainedMappedFixedValueFvPatchField
        //  onto a new patch
        selfContainedMappedFixedValueFvPatchField
        (
            const selfContainedMappedFixedValueFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        selfContainedMappedFixedValueFvPatchField
        (
            const selfContainedMappedFixedValueFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type> > clone() const
        {
            return tmp<fvPatchField<Type> >
            (
                new selfContainedMappedFixedValueFvPatchField<Type>
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        selfContainedMappedFixedValueFvPatchField
        (
            const selfContainedMappedFixedValueFvPatchField<Type>&,
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
                new selfContainedMappedFixedValueFvPatchField<Type>
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
#   include "selfContainedMappedFixedValueFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
