/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::codedFixedValueFvPatchField

Group
    grpGenericBoundaryConditions

Description
    Constructs on-the-fly a new boundary condition (derived from
    fixedValueFvPatchField) which is then used to evaluate.

    The code entries:
    \plaintable
       codeInclude | include files
       codeOptions | compiler line: added to EXE_INC (Make/options)
       codeLibs    | linker line: added to LIB_LIBS (Make/options)
       localCode   | c++; local static functions
       code        | c++; patch value assignment
       codeContext | additional dictionary context for the code
    \endplaintable

Usage
    Example:
    \verbatim
    <patchName>
    {
        type    codedFixedValue;
        value   uniform 0;
        name    rampedFixedValue;   // name of generated BC

        codeContext
        {
            ...
        }

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
    the code is read from a (runTimeModifiable!) dictionary system/codeDict
    which would have a corresponding entry:

    \verbatim
    <patchName>
    {
        code
        #{
            operator==(min(10, 0.1*this->db().time().value()));
        #};
    }
    \endverbatim

Note
    The code context dictionary can be supplied separately as the
    \c codeContext entry.

See also
    Foam::dynamicCode
    Foam::functionEntries::codeStream

SourceFiles
    codedFixedValueFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef codedFixedValueFvPatchField_H
#define codedFixedValueFvPatchField_H

#include "fixedValueFvPatchFields.H"
#include "codedBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class codedFixedValueFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class codedFixedValueFvPatchField
:
    public fixedValueFvPatchField<Type>,
    protected codedBase
{
    //- The parent boundary condition type
    typedef fixedValueFvPatchField<Type> parent_bctype;


    // Private Data

        //- Dictionary contents for the boundary condition
        dictionary dict_;

        const word name_;

        mutable autoPtr<fvPatchField<Type>> redirectPatchFieldPtr_;


    // Private Member Functions

        //- Mutable access to the loaded dynamic libraries
        virtual dlLibraryTable& libs() const;

        //- Description (type + name) for the output
        virtual string description() const;

        //- Clear redirected object(s)
        virtual void clearRedirect() const;

        //- Additional 'codeContext' dictionary to pass through
        virtual const dictionary& codeContext() const;

        //- The code dictionary. Inline "code" or from system/codeDict
        virtual const dictionary& codeDict() const;

        //- Adapt the context for the current object
        virtual void prepare(dynamicCode&, const dynamicCodeContext&) const;


public:

    // Static data members

        //- Name of the C code template to be used
        static constexpr const char* const codeTemplateC
            = "fixedValueFvPatchFieldTemplate.C";

        //- Name of the H code template to be used
        static constexpr const char* const codeTemplateH
            = "fixedValueFvPatchFieldTemplate.H";


    //- Runtime type information
    TypeName("codedFixedValue");


    // Constructors

        //- Construct from patch and internal field
        codedFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        codedFixedValueFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given codedFixedValueFvPatchField
        //  onto a new patch
        codedFixedValueFvPatchField
        (
            const codedFixedValueFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        codedFixedValueFvPatchField
        (
            const codedFixedValueFvPatchField<Type>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchField<Type>> clone() const
        {
            return tmp<fvPatchField<Type>>
            (
                new codedFixedValueFvPatchField<Type>(*this)
            );
        }

        //- Construct as copy setting internal field reference
        codedFixedValueFvPatchField
        (
            const codedFixedValueFvPatchField<Type>&,
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
                new codedFixedValueFvPatchField<Type>(*this, iF)
            );
        }


    // Member functions

        //- Get reference to the underlying patch
        const fvPatchField<Type>& redirectPatchField() const;

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Evaluate the patch field, sets Updated to false
        virtual void evaluate
        (
            const Pstream::commsTypes commsType=Pstream::commsTypes::blocking
        );

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "codedFixedValueFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
