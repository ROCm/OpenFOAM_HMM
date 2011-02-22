/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2011 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::codedFixedValueFvPatchScalarField

Description
    Constructs on-the-fly a new boundary condition (derived from
    fixedValueFvPatchScalarField) which is then used to evaluate.

    See also codeStream.

    Example:
    \verbatim
        movingWall
        {
            type            codedFixedValue<scalar>;
            value           uniform 0;
            redirectType    rampedFixedValue;   // name of generated bc

            code
            #{
                operator==(min(10, 0.1*this->db().time().value()));
            #};

            //codeInclude
            //#{
            //    #include "fvCFD.H"
            //#};

            //codeOptions
            //#{
            //    -I$(LIB_SRC)/finiteVolume/lnInclude
            //#};
        }
    \endverbatim

    A special form is if the 'code' section is not supplied. In this case
    the code gets read from a (runTimeModifiable!) dictionary system/codeDict
    which would have an entry

    \verbatim
    rampedFixedValue
    {
        code
        #{
            operator==(min(10, 0.1*this->db().time().value()));
        #};
    }
    \endverbatim

SeeAlso
    Foam::codeStreamTools for constant paths used

SourceFiles
    codedFixedValueFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef codedFixedValueFvPatchScalarField_H
#define codedFixedValueFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class codeProperties;

/*---------------------------------------------------------------------------*\
              Class codedFixedValueFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class codedFixedValueFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        mutable dictionary dict_;

        const word redirectType_;

        mutable autoPtr<fvPatchScalarField> redirectPatchFieldPtr_;


    // Private Member Functions

        const codeProperties& dict() const;

        void writeLibrary
        (
            const fileName& dir,
            const fileName& libPath,
            const dictionary& dict
        );

        void updateLibrary();

public:

    // Static data members

        //- Name of the C code template to be used
        const static word codeTemplateC;

        //- Name of the H code template to be used
        const static word codeTemplateH;


    //- Runtime type information
    TypeName("codedFixedValue<scalar>");


    // Constructors

        //- Construct from patch and internal field
        codedFixedValueFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        codedFixedValueFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  codedFixedValueFvPatchScalarField
        //  onto a new patch
        codedFixedValueFvPatchScalarField
        (
            const codedFixedValueFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        codedFixedValueFvPatchScalarField
        (
            const codedFixedValueFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new codedFixedValueFvPatchScalarField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        codedFixedValueFvPatchScalarField
        (
            const codedFixedValueFvPatchScalarField&,
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
                new codedFixedValueFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

        //- Get reference to the underlying patch
        const fvPatchScalarField& redirectPatchField() const;

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Evaluate the patch field, sets Updated to false
        virtual void evaluate
        (
            const Pstream::commsTypes commsType=Pstream::blocking
        );

        //- Write
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
