/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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
    Foam::pointPatchField

Description
    Abstract base class for point-mesh patch fields.

    The base-field does not store values as they are part of the
    "internal field".  There are derived classes to store constraint values
    e.g. fixedValuePointPatchField derived from the generic
    valuePointPatchField which ensures the values in the "internal field"
    are reset to the fixed-values by applying the stored values.

SourceFiles
    pointPatchField.C
    pointPatchFieldNew.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_pointPatchField_H
#define Foam_pointPatchField_H

#include "pointPatch.H"
#include "DimensionedField.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class objectRegistry;
class dictionary;
class pointPatchFieldMapper;
class pointMesh;

template<class Type> class pointPatchField;
template<class Type> class calculatedPointPatchField;

template<class Type>
Ostream& operator<<(Ostream&, const pointPatchField<Type>&);

/*---------------------------------------------------------------------------*\
                     Class pointPatchFieldBase Declaration
\*---------------------------------------------------------------------------*/

//- Template invariant parts for pointPatchField
class pointPatchFieldBase
{
    // Private Data

        //- Reference to patch
        const pointPatch& patch_;

        //- Update index used so that updateCoeffs is called only once during
        //- the construction of the matrix
        bool updated_;

        //- Optional patch type
        //  Used to allow specified boundary conditions to be applied
        //  to constraint patches by providing the constraint
        //  patch type as 'patchType'
        word patchType_;


protected:

    // Protected Member Functions

        //- Read dictionary entries.
        //  Useful when initially constructed without a dictionary
        virtual void readDict(const dictionary& dict);

        //- Set updated state
        void setUpdated(bool state) noexcept
        {
            updated_ = state;
        }


public:

    //- Debug switch to disallow the use of generic pointPatchField
    static int disallowGenericPatchField;

    //- Runtime type information
    TypeName("pointPatchField");


    // Constructors

        //- Construct from patch
        explicit pointPatchFieldBase(const pointPatch& p);

        //- Construct from patch and patch type
        pointPatchFieldBase(const pointPatch& p, const word& patchType);

        //- Construct from patch and dictionary
        pointPatchFieldBase(const pointPatch& p, const dictionary& dict);

        //- Copy construct with new patch
        pointPatchFieldBase(const pointPatchFieldBase&, const pointPatch& p);

        //- Copy construct
        pointPatchFieldBase(const pointPatchFieldBase&);


    //- Destructor
    virtual ~pointPatchFieldBase() = default;


    // Member Functions

    // Access

        //- The associated objectRegistry
        const objectRegistry& db() const;

        //- Return the patch
        const pointPatch& patch() const noexcept
        {
            return patch_;
        }

        //- The optional patch type
        const word& patchType() const noexcept
        {
            return patchType_;
        }

        //- The optional patch type
        word& patchType() noexcept
        {
            return patchType_;
        }


    // Solution

        //- True if the boundary condition has already been updated
        bool updated() const noexcept
        {
            return updated_;
        }


    // Check

        //- Check that patches are identical
        void checkPatch(const pointPatchFieldBase& rhs) const;
};


/*---------------------------------------------------------------------------*\
                       Class pointPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class pointPatchField
:
    public pointPatchFieldBase
{
    // Private Data

        //- Reference to internal field
        const DimensionedField<Type, pointMesh>& internalField_;


public:

    //- The Field value_type
    typedef Type value_type;

    //- The internal field type associated with the patch field
    typedef DimensionedField<Type, pointMesh> Internal;

    //- The patch type for the patch field
    typedef pointPatch Patch;

    //- Type for a \em calculated patch
    typedef calculatedPointPatchField<Type> Calculated;


    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            pointPatchField,
            patch,
            (
                const pointPatch& p,
                const DimensionedField<Type, pointMesh>& iF
            ),
            (p, iF)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            pointPatchField,
            patchMapper,
            (
                const pointPatchField<Type>& ptf,
                const pointPatch& p,
                const DimensionedField<Type, pointMesh>& iF,
                const pointPatchFieldMapper& m
            ),
            (dynamic_cast<const pointPatchFieldType&>(ptf), p, iF, m)
        );

        declareRunTimeSelectionTable
        (
            autoPtr,
            pointPatchField,
            dictionary,
            (
                const pointPatch& p,
                const DimensionedField<Type, pointMesh>& iF,
                const dictionary& dict
            ),
            (p, iF, dict)
        );


    // Constructors

        //- Construct from patch and internal field
        pointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Construct from patch, internal field and dictionary
        pointPatchField
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patchField<Type> onto a new patch
        pointPatchField
        (
            const pointPatchField<Type>&,
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Construct as copy
        pointPatchField(const pointPatchField<Type>&);

        //- Construct as copy setting internal field reference
        pointPatchField
        (
            const pointPatchField<Type>&,
            const DimensionedField<Type, pointMesh>&
        );

        //- Construct and return a clone
        virtual autoPtr<pointPatchField<Type>> clone() const = 0;

        //- Construct and return a clone setting internal field reference
        virtual autoPtr<pointPatchField<Type>> clone
        (
            const DimensionedField<Type, pointMesh>& iF
        ) const = 0;


    // Selectors

        //- Return a pointer to a new patchField created on freestore given
        //  patch and internal field
        //  (does not set the patch field values)
        static autoPtr<pointPatchField<Type>> New
        (
            const word& patchFieldType,
            const pointPatch& p,
            const DimensionedField<Type, pointMesh>& iF
        );

        //- Return a pointer to a new patchField created on freestore given
        //  patch and internal field
        //  (does not set the patch field values).
        //  Allows override of constraint type
        static autoPtr<pointPatchField<Type>> New
        (
            const word& patchFieldType,
            const word& actualPatchType,
            const pointPatch& p,
            const DimensionedField<Type, pointMesh>& iF
        );

        //- Return a pointer to a new patchField created on freestore from
        //  a given pointPatchField mapped onto a new patch
        static autoPtr<pointPatchField<Type>> New
        (
            const pointPatchField<Type>&,
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const pointPatchFieldMapper&
        );

        //- Return a pointer to a new patchField created on freestore
        //  from dictionary
        static autoPtr<pointPatchField<Type>> New
        (
            const pointPatch&,
            const DimensionedField<Type, pointMesh>&,
            const dictionary&
        );

        //- Return a pointer to a new calculatedPointPatchField created on
        //  freestore without setting patchField values
        template<class Type2>
        static autoPtr<pointPatchField<Type>>
        NewCalculatedType
        (
            const pointPatchField<Type2>&
        );


    //- Destructor
    virtual ~pointPatchField() = default;


    // Member Functions

        //- The type name for calculated patch fields
        static const word& calculatedType();

        //- The type name for zeroGradient patch fields
        static const word& zeroGradientType();


        // Attributes

            //- True if this patch field fixes a value
            virtual bool fixesValue() const
            {
                return false;
            }

            //- True if this patch field is coupled
            virtual bool coupled() const
            {
                return false;
            }

            //- The constraint type this pointPatchField implements.
            virtual const word& constraintType() const
            {
                return word::null;
            }


        // Access

            //- Return the patch size
            label size() const
            {
                return patch().size();
            }

            //- Return dimensioned internal field reference
            const DimensionedField<Type, pointMesh>& internalField()
            const noexcept
            {
                return internalField_;
            }

            //- Return internal field reference
            const Field<Type>& primitiveField() const noexcept
            {
                return internalField_;
            }

            //- Return field created from appropriate internal field values
            tmp<Field<Type>> patchInternalField() const;

            //- Return field created from appropriate internal field values
            //  given internal field reference
            template<class Type1>
            tmp<Field<Type1>> patchInternalField
            (
                const Field<Type1>& iF
            ) const;

            //- Return field created from selected internal field values
            //  given internal field reference
            template<class Type1>
            tmp<Field<Type1>> patchInternalField
            (
                const Field<Type1>& iF,
                const labelList& meshPoints
            ) const;

            //- Given the internal field and a patch field,
            //  add the patch field to the internal field
            template<class Type1>
            void addToInternalField
            (
                Field<Type1>& iF,
                const Field<Type1>& pF
            ) const;

            //- Given the internal field and a patch field,
            //  add selected elements of the patch field to the internal field
            template<class Type1>
            void addToInternalField
            (
                Field<Type1>& iF,
                const Field<Type1>& pF,
                const labelList& points
            ) const;

            //- Given the internal field and a patch field,
            //  set the patch field in the internal field
            template<class Type1>
            void setInInternalField
            (
                Field<Type1>& iF,
                const Field<Type1>& pF,
                const labelList& meshPoints
            ) const;

            //- Given the internal field and a patch field,
            //  set the patch field in the internal field
            template<class Type1>
            void setInInternalField
            (
                Field<Type1>& iF,
                const Field<Type1>& pF
            ) const;


        // Mapping Functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const pointPatchFieldMapper&
            )
            {}

            //- Reverse map the given pointPatchField onto this pointPatchField
            virtual void rmap
            (
                const pointPatchField<Type>&,
                const labelList&
            )
            {}


        // Evaluation Functions

            //- Update the coefficients associated with the patch field
            //  Sets Updated to true
            virtual void updateCoeffs();

            //- Initialise evaluation of the patch field (do nothing)
            virtual void initEvaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
            )
            {}

            //- Evaluate the patch field
            virtual void evaluate
            (
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::blocking
            );


        // I-O

            //- Write
            virtual void write(Ostream&) const;


    // Member Operators

        virtual void operator=(const pointPatchField<Type>&){}
        virtual void operator+=(const pointPatchField<Type>&){}
        virtual void operator-=(const pointPatchField<Type>&){}
        virtual void operator*=(const pointPatchField<scalar>&){}
        virtual void operator/=(const pointPatchField<scalar>&){}

        virtual void operator=(const Field<Type>&){}
        virtual void operator+=(const Field<Type>&){}
        virtual void operator-=(const Field<Type>&){}

        virtual void operator*=(const Field<scalar>&){}
        virtual void operator/=(const Field<scalar>&){}

        virtual void operator=(const Type&){}
        virtual void operator+=(const Type&){}
        virtual void operator-=(const Type&){}
        virtual void operator*=(const scalar){}
        virtual void operator/=(const scalar){}


        // Force an assignment irrespective of form of patch
        // By generic these do nothing unless the patch actually has boundary
        // values

            virtual void operator==(const pointPatchField<Type>&){}
            virtual void operator==(const Field<Type>&){}
            virtual void operator==(const Type&){}


    // Ostream Operator

        friend Ostream& operator<< <Type>
        (
            Ostream&,
            const pointPatchField<Type>&
        );


    // Other Methods

        //- Negate the field inplace. Dummy placeholder for FieldField
        void negate() {}

        //- Normalise the field inplace. Dummy placeholder for FieldField
        void normalise() {}
};


// This function is added as a hack to enable simple backward compatibility
// with versions using referenceLevel in GeometricField
template<class Type>
const pointPatchField<Type>& operator+
(
    const pointPatchField<Type>& ppf,
    const Type&
)
{
    return ppf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pointPatchFieldFunctions.H"

#ifdef NoRepository
    #include "pointPatchField.C"
    #include "calculatedPointPatchField.H"
    #include "zeroGradientPointPatchField.H"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Runtime selection macros

#define addToPointPatchFieldRunTimeSelection(PatchTypeField, typePatchTypeField) \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField,                                                    \
        patch                                                                  \
    );                                                                         \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField,                                                    \
        patchMapper                                                            \
    );                                                                         \
    addToRunTimeSelectionTable                                                 \
    (                                                                          \
        PatchTypeField,                                                        \
        typePatchTypeField,                                                    \
        dictionary                                                             \
    );


// for non-templated patch fields
#define makePointPatchTypeField(PatchTypeField,typePatchTypeField)             \
    defineTypeNameAndDebug(typePatchTypeField, 0);                             \
    addToPointPatchFieldRunTimeSelection(PatchTypeField, typePatchTypeField)


// for templated patch fields
#define makeTemplatePointPatchTypeField(PatchTypeField, typePatchTypeField)    \
    defineNamedTemplateTypeNameAndDebug(typePatchTypeField, 0);                \
    addToPointPatchFieldRunTimeSelection(PatchTypeField, typePatchTypeField)


#define makePointPatchFields(type)                                             \
    makeTemplatePointPatchTypeField                                            \
    (                                                                          \
        pointPatchScalarField,                                                 \
        type##PointPatchScalarField                                            \
    );                                                                         \
    makeTemplatePointPatchTypeField                                            \
    (                                                                          \
        pointPatchVectorField,                                                 \
        type##PointPatchVectorField                                            \
    );                                                                         \
    makeTemplatePointPatchTypeField                                            \
    (                                                                          \
        pointPatchSphericalTensorField,                                        \
        type##PointPatchSphericalTensorField                                   \
    );                                                                         \
    makeTemplatePointPatchTypeField                                            \
    (                                                                          \
        pointPatchSymmTensorField,                                             \
        type##PointPatchSymmTensorField                                        \
    );                                                                         \
    makeTemplatePointPatchTypeField                                            \
    (                                                                          \
        pointPatchTensorField,                                                 \
        type##PointPatchTensorField                                            \
    );


#define makePointPatchFieldsTypeName(type)                                     \
    defineNamedTemplateTypeNameAndDebug(type##PointPatchScalarField, 0);       \
    defineNamedTemplateTypeNameAndDebug(type##PointPatchVectorField, 0);       \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        type##PointPatchSphericalTensorField, 0                                \
    );                                                                         \
    defineNamedTemplateTypeNameAndDebug(type##PointPatchSymmTensorField, 0);   \
    defineNamedTemplateTypeNameAndDebug(type##PointPatchTensorField, 0)


#define makePointPatchFieldTypedefs(type)                                      \
    typedef type##PointPatchField<scalar> type##PointPatchScalarField;         \
    typedef type##PointPatchField<vector> type##PointPatchVectorField;         \
    typedef type##PointPatchField<sphericalTensor>                             \
        type##PointPatchSphericalTensorField;                                  \
    typedef type##PointPatchField<symmTensor> type##PointPatchSymmTensorField;\
    typedef type##PointPatchField<tensor> type##PointPatchTensorField;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
