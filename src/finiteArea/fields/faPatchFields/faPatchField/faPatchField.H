/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::faPatchField

Description
    faPatchField<Type> abstract base class.  This class gives a fat-interface
    to all derived classes covering all possible ways in which they might be
    used.  The first level of derivation is to basic patchFields which cover
    zero-gradient, fixed-gradient, fixed-value and mixed conditions.  The next
    level of derivation covers all the specialised typed with specific
    evaluation procedures, particularly with respect to specific fields.

Author
    Zeljko Tukovic, FMENA
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    faPatchField.C
    newPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef faPatchField_H
#define faPatchField_H

#include "faPatch.H"
#include "DimensionedField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations

class objectRegistry;
class dictionary;
class faPatchFieldMapper;
class areaMesh;

template<class Type>
class faPatchField;

template<class Type>
class calculatedFaPatchField;

template<class Type>
Ostream& operator<<(Ostream&, const faPatchField<Type>&);

/*---------------------------------------------------------------------------*\
                        Class faPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class faPatchField
:
    public Field<Type>
{
    // Private Data

        //- Reference to a patch
        const faPatch& patch_;

        //- Reference to internal field
        const DimensionedField<Type, areaMesh>& internalField_;

        //- Update index used so that updateCoeffs is called only once during
        //  the construction of the matrix
        bool updated_;

        //- Optional patch type, used to allow specified boundary conditions
        //- to be applied to constraint patches by providing the constraint
        //- patch type as 'patchType'
        word patchType_;


public:

    typedef faPatch Patch;
    typedef calculatedFaPatchField<Type> Calculated;


    //- Runtime type information
    TypeName("faPatchField");

    //- Debug switch to disallow the use of
    static int disallowGenericFaPatchField;


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            tmp,
            faPatchField,
            patch,
            (
                const faPatch& p,
                const DimensionedField<Type, areaMesh>& iF
            ),
            (p, iF)
        );

        declareRunTimeSelectionTable
        (
            tmp,
            faPatchField,
            patchMapper,
            (
                const faPatchField<Type>& ptf,
                const faPatch& p,
                const DimensionedField<Type, areaMesh>& iF,
                const faPatchFieldMapper& m
            ),
            (dynamic_cast<const faPatchFieldType&>(ptf), p, iF, m)
        );

        declareRunTimeSelectionTable
        (
            tmp,
            faPatchField,
            dictionary,
            (
                const faPatch& p,
                const DimensionedField<Type, areaMesh>& iF,
                const dictionary& dict
            ),
            (p, iF, dict)
        );


    // Constructors

        //- Construct from patch and internal field
        faPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&
        );

        //- Construct from patch and internal field and patch field
        faPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const Field<Type>&
        );

        //- Construct from patch, internal field and dictionary
        faPatchField
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const dictionary&
        );

        //- Construct by mapping the given faPatchField onto a new patch
        faPatchField
        (
            const faPatchField<Type>&,
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const faPatchFieldMapper&
        );

        //- Construct as copy
        faPatchField(const faPatchField<Type>&);

        //- Construct and return a clone
        virtual tmp<faPatchField<Type>> clone() const
        {
            return tmp<faPatchField<Type>>(new faPatchField<Type>(*this));
        }

        //- Construct as copy setting internal field reference
        faPatchField
        (
            const faPatchField<Type>&,
            const DimensionedField<Type, areaMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<faPatchField<Type>> clone
        (
            const DimensionedField<Type, areaMesh>& iF
        ) const
        {
            return tmp<faPatchField<Type>>(new faPatchField<Type>(*this, iF));
        }


    // Selectors

        //- Return a pointer to a new patchField created on freestore given
        //- patch and internal field
        //  (does not set the patch field values)
        static tmp<faPatchField<Type>> New
        (
            const word& patchFieldType,
            const word& actualPatchType,
            const faPatch&,
            const DimensionedField<Type, areaMesh>&
        );

        //- Return a pointer to a new patchField created on freestore given
        //- patch and internal field
        //  (does not set the patch field values)
        static tmp<faPatchField<Type>> New
        (
            const word& patchFieldType,
            const faPatch&,
            const DimensionedField<Type, areaMesh>&
        );

        //- Return a pointer to a new patchField created on freestore from
        //- a given faPatchField mapped onto a new patch
        static tmp<faPatchField<Type>> New
        (
            const faPatchField<Type>&,
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const faPatchFieldMapper&
        );

        //- Return a pointer to a new patchField created on freestore
        //- from dictionary
        static tmp<faPatchField<Type>> New
        (
            const faPatch&,
            const DimensionedField<Type, areaMesh>&,
            const dictionary&
        );

        //- Return a pointer to a new calculatedFaPatchField created on
        //- freestore without setting patchField values
        template<class Type2>
        static tmp<faPatchField<Type>> NewCalculatedType
        (
            const faPatchField<Type2>&
        );


    //- Destructor
    virtual ~faPatchField<Type>() = default;


    // Member Functions

        // Access

            //- Return local objectRegistry
            const objectRegistry& db() const;

            //- Return patch
            const faPatch& patch() const
            {
                return patch_;
            }

            //- Return dimensioned internal field reference
            const DimensionedField<Type, areaMesh>& internalField() const
            {
                return internalField_;
            }

            //- Return internal field reference
            const Field<Type>& primitiveField() const
            {
                return internalField_;
            }

            //- Optional patch type
            const word& patchType() const
            {
                return patchType_;
            }

            //- Optional patch type
            word& patchType()
            {
                return patchType_;
            }

            //- Return the type of the calculated for of faPatchField
            static const word& calculatedType();

            //- Return true if this patch field fixes a value.
            //  Needed to check if a level has to be specified while solving
            //  Poissons equations.
            virtual bool fixesValue() const
            {
                return false;
            }

            //- Return true if this patch field is coupled
            virtual bool coupled() const
            {
                return false;
            }

            //- Return true if the boundary condition has already been updated
            bool updated() const
            {
                return updated_;
            }


        // Mapping

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const faPatchFieldMapper&
            );

            //- Reverse map the given faPatchField onto this faPatchField
            virtual void rmap
            (
                const faPatchField<Type>&,
                const labelList&
            );


        // Evaluation

            //- Return patch-normal gradient
            virtual tmp<Field<Type>> snGrad() const;

            //- Update the coefficients associated with the patch field
            //  Sets Updated to true
            virtual void updateCoeffs()
            {
                updated_ = true;
            }

            //- Return internal field next to patch as patch field
            virtual tmp<Field<Type>> patchInternalField() const;

            //- Return patchField on the opposite patch of a coupled patch
            virtual tmp<Field<Type>> patchNeighbourField() const
            {
                NotImplemented;
                return *this;
            }

            //- Initialise the evaluation of the patch field
            virtual void initEvaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
            )
            {}

            //- Evaluate the patch field, sets Updated to false
            virtual void evaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
            );

            //- Return the matrix diagonal coefficients corresponding to the
            //- evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type>> valueInternalCoeffs
            (
                const tmp<Field<scalar>>&
            ) const
            {
                NotImplemented;
                return *this;
            }

            //- Return the matrix source coefficients corresponding to the
            //- evaluation of the value of this patchField with given weights
            virtual tmp<Field<Type>> valueBoundaryCoeffs
            (
                const tmp<Field<scalar>>&
            ) const
            {
                NotImplemented;
                return *this;
            }

            //- Return the matrix diagonal coefficients corresponding to the
            //- evaluation of the gradient of this patchField
            virtual tmp<Field<Type>> gradientInternalCoeffs() const
            {
                NotImplemented;
                return *this;
            }

            //- Return the matrix source coefficients corresponding to the
            //- evaluation of the gradient of this patchField
            virtual tmp<Field<Type>> gradientBoundaryCoeffs() const
            {
                NotImplemented;
                return *this;
            }


        // IO

            //- Write
            virtual void write(Ostream&) const;


        // Check

            //- Check faPatchField<Type> against given faPatchField<Type>
            void check(const faPatchField<Type>&) const;


    // Member Operators

        virtual void operator=(const UList<Type>&);

        virtual void operator=(const faPatchField<Type>&);
        virtual void operator+=(const faPatchField<Type>&);
        virtual void operator-=(const faPatchField<Type>&);
        virtual void operator*=(const faPatchField<scalar>&);
        virtual void operator/=(const faPatchField<scalar>&);

        virtual void operator+=(const Field<Type>&);
        virtual void operator-=(const Field<Type>&);

        virtual void operator*=(const Field<scalar>&);
        virtual void operator/=(const Field<scalar>&);

        virtual void operator=(const Type&);
        virtual void operator+=(const Type&);
        virtual void operator-=(const Type&);
        virtual void operator*=(const scalar);
        virtual void operator/=(const scalar);


        // Force an assignment irrespective of form of patch

        virtual void operator==(const faPatchField<Type>&);
        virtual void operator==(const Field<Type>&);
        virtual void operator==(const Type&);


    // Ostream operator

        friend Ostream& operator<< <Type>(Ostream&, const faPatchField<Type>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "faPatchField.C"
    #include "calculatedFaPatchField.H"
#endif


#define addToFaPatchFieldRunTimeSelection(PatchTypeField, typePatchTypeField) \
                                                                              \
addToRunTimeSelectionTable                                                    \
(                                                                             \
    PatchTypeField, typePatchTypeField, patch                                 \
);                                                                            \
                                                                              \
addToRunTimeSelectionTable                                                    \
(                                                                             \
    PatchTypeField,                                                           \
    typePatchTypeField,                                                       \
    patchMapper                                                               \
);                                                                            \
                                                                              \
addToRunTimeSelectionTable                                                    \
(                                                                             \
    PatchTypeField, typePatchTypeField, dictionary                            \
);


#define makeFaPatchTypeFieldTypeName(typePatchTypeField)                      \
                                                                              \
defineNamedTemplateTypeNameAndDebug(typePatchTypeField, 0);


#define makeFaPatchFieldsTypeName(typePatchField)                             \
                                                                              \
makeFaPatchTypeFieldTypeName(typePatchField##FaPatchScalarField);             \
makeFaPatchTypeFieldTypeName(typePatchField##FaPatchVectorField);             \
makeFaPatchTypeFieldTypeName(typePatchField##FaPatchSphericalTensorField);    \
makeFaPatchTypeFieldTypeName(typePatchField##FaPatchSymmTensorField);         \
makeFaPatchTypeFieldTypeName(typePatchField##FaPatchTensorField);


#define makeFaPatchTypeField(PatchTypeField, typePatchTypeField)              \
                                                                              \
defineTypeNameAndDebug(typePatchTypeField, 0);                                \
                                                                              \
addToFaPatchFieldRunTimeSelection                                             \
(                                                                             \
    PatchTypeField, typePatchTypeField                                        \
);

#define makeTemplateFaPatchTypeField(PatchTypeField, typePatchTypeField)      \
                                                                              \
defineNamedTemplateTypeNameAndDebug(typePatchTypeField, 0);                   \
                                                                              \
addToFaPatchFieldRunTimeSelection                                             \
(                                                                             \
    PatchTypeField, typePatchTypeField                                        \
);


#define makeFaPatchFields(type)                                               \
                                                                              \
makeTemplateFaPatchTypeField(faPatchScalarField, type##FaPatchScalarField);   \
makeTemplateFaPatchTypeField(faPatchVectorField, type##FaPatchVectorField);   \
makeTemplateFaPatchTypeField                                                  \
(                                                                             \
    faPatchSphericalTensorField,                                              \
    type##FaPatchSphericalTensorField                                         \
);                                                                            \
makeTemplateFaPatchTypeField                                                  \
(                                                                             \
    faPatchSymmTensorField,                                                   \
    type##FaPatchSymmTensorField                                              \
);                                                                            \
makeTemplateFaPatchTypeField                                                  \
(                                                                             \
    faPatchTensorField,                                                       \
    type##FaPatchTensorField                                                  \
);


#define makeFaPatchTypeFieldTypedefs(type)                                    \
                                                                              \
typedef type##FaPatchField<scalar> type##FaPatchScalarField;                  \
typedef type##FaPatchField<vector> type##FaPatchVectorField;                  \
typedef type##FaPatchField<sphericalTensor>                                   \
     type##FaPatchSphericalTensorField;                                       \
typedef type##FaPatchField<symmTensor> type##FaPatchSymmTensorField;          \
typedef type##FaPatchField<tensor> type##FaPatchTensorField;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
