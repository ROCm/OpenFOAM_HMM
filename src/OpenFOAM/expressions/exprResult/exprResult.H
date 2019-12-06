/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2018 Bernhard Gschaider <bgschaid@hfd-research.com>
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
    Foam::expressions::exprResult

Description
    A polymorphic field/result from evaluating an expression

    \heading Dictionary parameters
    \table
        Property    | Description                           | Required | Default
        resultType  | The type of result                    | no  | exprResult
        unsetValue  | Create without reading the dictionary | no  | false
        noReset     | Suppress reset on time                | no  | false
    \endtable

    When creating with values
    \table
        Property    | Description                           | Required | Default
        valueType   | Result value type (scalar, vector,..) | yes |
        isSingleValue | A single field value                | no  | false
        isPointValue  | Interpret values as point values    | no  | false
        value       | The field values                      | yes |
        fieldSize   | The size of the field (when not single-value) | no  |
    \endtable

SourceFiles
    exprResult.C
    exprResultI.H

\*---------------------------------------------------------------------------*/

#ifndef expressions_exprResult_H
#define expressions_exprResult_H

#include "vector.H"
#include "tensor.H"
#include "sphericalTensor.H"
#include "symmTensor.H"
#include "dimensionedType.H"
#include "IOField.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace expressions
{

/*---------------------------------------------------------------------------*\
                      Class exprResult Declaration
\*---------------------------------------------------------------------------*/

class exprResult
:
    public refCount
{
    // Private Data

        //- The value type as string,
        //- normally corresponds to pTraits or typeName
        word valType_;

        //- Is single, uniform value (can be a non-field)
        bool isUniform_;

        //- Is a point value
        bool isPointVal_;

        //- Whether or not the variable will be reset
        bool noReset_;

        //- Allow override of noReset_, but only accessible for subclasses
        bool needsReset_;

        //- Size of field or object
        label size_;

        //- A %union of single values, including standard VectorSpace types
        union singleValue
        {
            bool bool_;
            label label_;
            scalar scalar_;
            vector vector_;
            tensor tensor_;
            symmTensor symmTensor_;
            sphericalTensor sphTensor_;

            //- Construct null, zero-initialized
            singleValue();

            //- Copy construct
            singleValue(const singleValue& val);

            //- Copy assignment
            void operator=(const singleValue& val);

            //- Return current value for specified type.
            template<class T>
            inline const T& get() const
            {
                WarningInFunction
                    << "Not implemented for type " << pTraits<T>::typeName << nl;
                return pTraits<T>::zero;
            }

            //- Set new value for specified type. Returns updated value.
            template<class T>
            inline const T& set(const T& val)
            {
                WarningInFunction
                    << "Not implemented for type " << pTraits<T>::typeName << nl;
                return val;
            }
        };

        //- The single value
        singleValue single_;

        //- Allocated plain field (eg, scalarField)
        void *fieldPtr_;

        //- Alternative storage for non-plain fields (eg, volScalarField)
        autoPtr<regIOobject> objectPtr_;


    // Private Member Functions

        //- Type-checked deletion of the value pointer.
        //  \return True if the type check was satisfied
        template<class Type>
        inline bool deleteChecked();

        //- Dispatch to type-checked pointer deletion
        void uglyDelete();

        //- Type-checked creation of field from dictionary
        //  \return True if the type check was satisfied
        template<class Type>
        inline bool readChecked
        (
            const word& key,
            const dictionary& dict,
            const label len,
            const bool uniform
        );

        //- Type-checked retrieval of uniform field from current results
        //  \return True if the type check was satisfied
        template<class Type>
        bool getUniformChecked
        (
            exprResult& result,
            const label size,
            const bool noWarn,
            const bool parRun
        ) const;

        //- Type-checked retrieval of uniform field from current results
        //  \return True if the type check was satisfied
        bool getUniformCheckedBool
        (
            exprResult& result,
            const label size,
            const bool noWarn,
            const bool parRun
        ) const;

        //- Type-checked determination of centre value (min/max)
        //  \return True if the type check was satisfied
        template<class Type>
        bool setAverageValueChecked(const bool parRun = Pstream::parRun());

        //- Type-checked determination of average bool value
        //  \return True if the type check was satisfied
        bool setAverageValueCheckedBool(const bool parRun = Pstream::parRun());

        //- Type-checked copy of field
        //  \return True if the type check was satisfied
        template<class Type>
        bool duplicateFieldChecked(const void* ptr);

        //- Type-checked writing of the single value (uniform) entry
        //  \return True if the type check was satisfied
        template<class Type>
        bool writeSingleValueChecked(Ostream& os) const;

        //- Type-checked writing of "value" field entry
        //  \return True if the type check was satisfied
        template<class Type>
        bool writeValueFieldChecked(Ostream& os) const;

        //- Type-checked forwarding to Field::writeEntry
        //  \return True if the type check was satisfied
        template<class Type>
        bool writeEntryChecked(const word& keyword, Ostream& os) const;

        //- Type-checked field addition with another expression field
        //  \return True if the type check was satisfied
        template<class Type>
        bool plusEqChecked(const exprResult& b);

        //- Type-checked field multiplication with a scalar
        //  \return True if the type check was satisfied
        template<class Type>
        bool multiplyEqChecked(const scalar& b);


        template<class Type>
        inline void setResultImpl(Field<Type>*, bool isPointVal=false);

        template<class Type>
        inline void setResultImpl(const Field<Type>&, bool isPointVal=false);

        template<class Type>
        inline void setResultImpl(Field<Type>&&, bool isPointVal=false);

        template<class Type>
        inline void setResultImpl(const Type& val, const label len);

        template<class Type>
        inline void setSingleValueImpl(const Type& val);

        template<class Type>
        inline void setObjectResultImpl(Type* ptr);

        template<class Type>
        inline void setObjectResultImpl(autoPtr<Type>& ap);

        template<class Type>
        inline void setObjectResultImpl(autoPtr<Type>&& ap);


protected:

    // Protected Member Functions

        //- Simulate virtual templated methods
        inline virtual exprResult& target() { return *this; }

        //- Reset at new timestep according to the derived class type
        virtual void resetImpl();

        //- Reset at new timestep according to type
        //  \return true if it was actually reset
        bool reset(bool force=false);

        //- Adjusts the internal needsReset value
        void needsReset(bool val) { needsReset_ = val; }


public:

        //- An empty result
        static const exprResult null;

        //- Friendship with globals
        friend class exprResultGlobals;


        //- Runtime type information
        TypeName("exprResult");

        declareRunTimeSelectionTable
        (
            autoPtr,
            exprResult,
            dictionary,
            (
                const dictionary& dict
            ),
            (dict)
        );
        declareRunTimeSelectionTable
        (
            autoPtr,
            exprResult,
            empty,
            (),
            ()
        );


    // Constructors

        //- Construct null
        exprResult();

        //- Copy construct
        exprResult(const exprResult& expr);

        //- Move construct
        exprResult(exprResult&& expr);

        //- Construct from a dictionary
        explicit exprResult
        (
            const dictionary& dict,
            const bool uniform = false,
            const bool needsValue = false
        );

        //- Construct by copying a field
        template<class Type>
        exprResult(const Field<Type>& f);

        //- Construct by moving a field
        template<class Type>
        exprResult(Field<Type>&& f);

        //- Construct for an IOobject
        template<class Type>
        exprResult(autoPtr<Type>& ap);

        //- Construct for an IOobject
        template<class Type>
        exprResult(autoPtr<Type>&& ap);

        //- Construct from a dimensioned value
        template<class Type>
        exprResult(const dimensioned<Type>& f);

        #undef exprResult_Construct
        #define exprResult_Construct(Type)                                \
            /*! Construct from single value of Type */                    \
            explicit exprResult(const Type& val) : exprResult()           \
            {                                                             \
                setSingleValue(val);                                      \
            }

        exprResult_Construct(bool);
        exprResult_Construct(scalar);
        exprResult_Construct(vector);
        exprResult_Construct(tensor);
        exprResult_Construct(symmTensor);
        exprResult_Construct(sphericalTensor);

        #undef exprResult_Construct


    // Selectors

        //- Return a reference to the selected value driver
        static autoPtr<exprResult> New(const dictionary& dict);

        //- Clone
        virtual autoPtr<exprResult> clone() const
        {
            return autoPtr<exprResult>::New(*this);
        }


    //- Destructor
    virtual ~exprResult();


    // Member Functions

    // Access

        //- Has a value?
        inline bool hasValue() const;

        //- Basic type for the field or single value
        inline const word& valueType() const;

        //- True if representing point values, or test if same as isPointVal
        inline bool isPointValue(const bool isPointVal = true) const;

        //- True if single, uniform value
        inline bool isUniform() const;

        //- True if valueType corresponds to the given Type
        template<class Type>
        inline bool isType() const;

        //- True if valueType is a bool
        inline bool isBool() const;

        //- True if the object pointer is being used
        inline bool isObject() const;

        //- The field or object size
        inline label size() const;

        //- The address of the field data content.
        //  Fatal for unknown types.
        //  Used, for example, for python integration
        const void* dataAddress() const;


    // Edit

        //- Clear (zero) the result
        void clear();

        //- Change reset behaviour
        void noReset() { noReset_ = true; }

        //- Change reset behaviour
        void allowReset() { noReset_ = false; }

        //- Test if field corresponds to a single-value and thus uniform.
        //  Uses field min/max to establish uniformity.
        //  Test afterwards with isUniform()
        void testIfSingleValue(const bool parRun = Pstream::parRun());


    // Set results

        //- Set result field, taking ownership of the pointer
        template<class Type>
        inline void setResult(Field<Type>*, bool isPointVal=false);

        //- Set result field, taking copy of the field contents
        template<class Type>
        inline void setResult(const Field<Type>& fld, bool isPointVal=false);

        //- Set result field, moving field contents
        template<class Type>
        inline void setResult(Field<Type>&&, bool isPointVal=false);

        //- Set uniform result field of given size
        template<class Type>
        inline void setResult(const Type& val, const label size);

        //- Set single-value uniform result
        template<class Type>
        inline void setSingleValue(const Type& val);

        template<class Type>
        inline void setObjectResult(autoPtr<Type>& o);

        template<class Type>
        inline void setObjectResult(autoPtr<Type>&& o);


    // Access/Get results

        //- Return const reference to the field
        template<class Type>
        inline const Field<Type>& cref() const;

        //- Return non-const reference to the field
        template<class Type>
        inline Field<Type>& ref();

        //- Return non-const reference to the field, casting away constness
        template<class Type>
        inline Field<Type>& getRef() const;

        //- Return tmp field of the contents,
        //- optionally keeping a copy in cache
        template<class Type>
        inline tmp<Field<Type>> getResult(bool cacheCopy=false);

        //- Get object result (Caution - potentially fragile)
        //- optionally keeping a copy in cache
        template<class Type>
        inline tmp<Type> getObjectResult(bool cacheCopy=false);

        //- Construct a uniform field from the current results
        //  Uses the field average. Optionally warning if the min/max
        //  deviation is larger than SMALL.
        exprResult getUniform
        (
            const label size,
            const bool noWarn,
            const bool parRun = Pstream::parRun()
        ) const;

        //- Get a reduced result
        template<template<class> class BinaryOp, class Type>
        inline Type getReduced
        (
            const BinaryOp<Type>& bop,
            const Type& initial = pTraits<Type>::zero
        );


    // Write

        //- Forwarding to Field::writeEntry
        void writeEntry(const word& keyword, Ostream& os) const;

        //- Write entry as dictionary contents
        void writeDict(Ostream& os, const bool subDict=true) const;

        //- Write the single value, or the first value from field
        void writeValue(Ostream& os) const;


    // Member Operators

        //- Copy assignment
        virtual void operator=(const exprResult& rhs);

        //- Move assignment
        virtual void operator=(exprResult&& rhs);


        //- Scalar multiplication
        exprResult& operator*=(const scalar& b);

        //- Addition of results
        exprResult& operator+=(const exprResult& b);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace expressions

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Operators

expressions::exprResult operator*
(
    const scalar& a,
    const expressions::exprResult& b
);
expressions::exprResult operator*
(
    const expressions::exprResult& a,
    const scalar& b
);
expressions::exprResult operator+
(
    const expressions::exprResult& a,
    const expressions::exprResult& b
);


// IO Operator
Istream& operator>>(Istream& os, expressions::exprResult& data);
Ostream& operator<<(Ostream& os, const expressions::exprResult& data);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "exprResultI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
