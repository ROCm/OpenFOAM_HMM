/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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
    Foam::GeometricField

Description
    Generic GeometricField class.

SourceFiles
    GeometricFieldI.H
    GeometricField.C
    GeometricBoundaryField.C
    GeometricFieldFunctions.H
    GeometricFieldFunctions.C

\*---------------------------------------------------------------------------*/

#ifndef GeometricField_H
#define GeometricField_H

#include "regIOobject.H"
#include "dimensionedTypes.H"
#include "DimensionedField.H"
#include "FieldField.H"
#include "lduInterfaceFieldPtrsList.H"
#include "LduInterfaceFieldPtrsList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class dictionary;

template<class Type, template<class> class PatchField, class GeoMesh>
class GeometricField;

template<class Type, template<class> class PatchField, class GeoMesh>
Ostream& operator<<
(
    Ostream&,
    const GeometricField<Type, PatchField, GeoMesh>&
);

template<class Type, template<class> class PatchField, class GeoMesh>
Ostream& operator<<
(
    Ostream&,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>&
);


/*---------------------------------------------------------------------------*\
                           Class GeometricField Declaration
\*---------------------------------------------------------------------------*/

template<class Type, template<class> class PatchField, class GeoMesh>
class GeometricField
:
    public DimensionedField<Type, GeoMesh>
{
    // Private Member Functions

        //- Read from file if it is present
        bool readIfPresent();

        //- Read old time field from file if it is present
        bool readOldTimeIfPresent();


public:

    // Public Typedefs

        //- Type of mesh on which this GeometricField is instantiated
        typedef typename GeoMesh::Mesh Mesh;

        //- Type of boundary mesh on which this GeometricField::Boundary is
        //- instantiated
        typedef typename GeoMesh::BoundaryMesh BoundaryMesh;

        //- Type of the internal field from which this GeometricField is derived
        typedef DimensionedField<Type, GeoMesh> Internal;

        //- Type of the patch field of which the GeometricField::Boundary is
        //- composed
        typedef PatchField<Type> Patch;


    //- The boundary fields
    class Boundary
    :
        public FieldField<PatchField, Type>
    {
        // Private Data

            //- Reference to BoundaryMesh for which this field is defined
            const BoundaryMesh& bmesh_;


    public:

        // Constructors

            //- Construct from a BoundaryMesh, setting patches later
            explicit Boundary(const BoundaryMesh& bmesh);

            //- Construct from a BoundaryMesh, reference to the internal field
            //- and a patch type
            Boundary
            (
                const BoundaryMesh& bmesh,
                const Internal& field,
                const word& patchFieldType = PatchField<Type>::calculatedType()
            );

            //- Construct from a BoundaryMesh, reference to the internal field
            //- and a wordList of patch types and optional the actual patch
            //- types (to override constraint patches)
            Boundary
            (
                const BoundaryMesh& bmesh,
                const Internal& field,
                const wordList& wantedPatchTypes,
                const wordList& actualPatchTypes = wordList()
            );

            //- Construct from a BoundaryMesh, reference to the internal field
            //- and a PtrList<PatchField<Type>>
            Boundary
            (
                const BoundaryMesh& bmesh,
                const Internal& field,
                const PtrList<PatchField<Type>>&
            );

            //- Construct as copy setting the reference to the internal field
            Boundary
            (
                const Internal& field,
                const Boundary& btf
            );

            //- Construct as copy setting the reference to the internal field
            //- and resetting type of field for given patch IDs
            Boundary
            (
                const Internal& field,
                const Boundary& btf,
                const labelList& patchIDs,
                const word& patchFieldName
            );

            //- Copy construct
            //  Dangerous because Field may be set to a field which gets deleted
            //  Need new type of BoundaryField, one which is part of a geometric
            //  field for which snGrad etc. may be called and a free standing
            //  BoundaryField for which such operations are unavailable.
            Boundary(const Boundary& btf);

            //- Construct from dictionary
            Boundary
            (
                const BoundaryMesh& bmesh,
                const Internal& field,
                const dictionary& dict
            );


        // Member Functions

            //- Read the boundary field
            void readField
            (
                const Internal& field,
                const dictionary& dict
            );

            //- Update the boundary condition coefficients
            void updateCoeffs();

            //- Evaluate boundary conditions
            void evaluate();

            //- Return a list of the patch types
            wordList types() const;

            //- Return BoundaryField of the cell values neighbouring
            //- the boundary
            Boundary boundaryInternalField() const;

            //- Return a list of pointers for each patch field with only those
            //- pointing to interfaces being set
            LduInterfaceFieldPtrsList<Type> interfaces() const;

            //- Return a list of pointers for each patch field with only those
            //- pointing to interfaces being set
            lduInterfaceFieldPtrsList scalarInterfaces() const;

            //- Write boundary field as dictionary entry
            void writeEntry(const word& keyword, Ostream& os) const;

            //- Write dictionary entries of the individual boundary fields.
            void writeEntries(Ostream& os) const;


        // Member Operators

            //- Assignment to BoundaryField<Type, PatchField, BoundaryMesh>
            void operator=(const Boundary&);

            //- Assignment to FieldField<PatchField, Type>
            void operator=(const FieldField<PatchField, Type>&);

            //- Assignment to Type
            void operator=(const Type&);


            //- Forced assignment to
            //  BoundaryField<Type, PatchField, BoundaryMesh>
            void operator==(const Boundary&);

            //- Forced assignment to FieldField<PatchField, Type>
            void operator==(const FieldField<PatchField, Type>&);

            //- Forced assignment to Type
            void operator==(const Type&);
    };


private:

    // Private data

        //- Current time index.
        //  Used to trigger the storing of the old-time value
        mutable label timeIndex_;

        //- Pointer to old time field
        mutable GeometricField<Type, PatchField, GeoMesh>* field0Ptr_;

        //- Pointer to previous iteration (used for under-relaxation)
        mutable GeometricField<Type, PatchField, GeoMesh>* fieldPrevIterPtr_;

        //- Boundary Type field containing boundary field values
        Boundary boundaryField_;


    // Private Member Functions

        //- Read the field from the dictionary
        void readFields(const dictionary& dict);

        //- Read the field - create the field dictionary on-the-fly
        void readFields();


public:

    //- Runtime type information
    TypeName("GeometricField");


    // Public typedefs

        typedef typename Field<Type>::cmptType cmptType;

    // Static Member Functions

        //- Return a null geometric field
        inline static const GeometricField<Type, PatchField, GeoMesh>& null();


    // Constructors

        //- Construct given IOobject, mesh, dimensions and patch type.
        //  This allocates storage for the field but does not set values.
        //  Used only within this class to create TEMPORARY variables
        GeometricField
        (
            const IOobject& io,
            const Mesh& mesh,
            const dimensionSet& ds,
            const word& patchFieldType = PatchField<Type>::calculatedType()
        );

        //- Construct given IOobject, mesh, dimensions and patch types.
        //  This allocates storage for the field but does not set values.
        //  Used only within this class to create TEMPORARY variables
        GeometricField
        (
            const IOobject& io,
            const Mesh& mesh,
            const dimensionSet& ds,
            const wordList& wantedPatchTypes,
            const wordList& actualPatchTypes = wordList()
        );

        //- Construct given IOobject, mesh, dimensioned<Type> and patch type.
        //  This assigns both dimensions and values.
        //  The internal name for the dimensioned\<Type\> has no influence.
        GeometricField
        (
            const IOobject& io,
            const Mesh& mesh,
            const dimensioned<Type>& dt,
            const word& patchFieldType = PatchField<Type>::calculatedType()
        );

        //- Construct given IOobject, mesh, dimensioned<Type> and patch types.
        //  This assigns both dimensions and values.
        //  The internal name for the dimensioned\<Type\> has no influence.
        GeometricField
        (
            const IOobject& io,
            const Mesh& mesh,
            const dimensioned<Type>& dt,
            const wordList& wantedPatchTypes,
            const wordList& actualPatchTypes = wordList()
        );

        //- Copy construct from components
        GeometricField
        (
            const IOobject& io,
            const Internal& diField,
            const PtrList<PatchField<Type>>& ptfl
        );

        //- Copy construct from internal field, with specified patch type
        GeometricField
        (
            const IOobject& io,
            const Mesh& mesh,
            const dimensionSet& ds,
            const Field<Type>& iField,
            const word& patchFieldType = PatchField<Type>::calculatedType()
        );

        //- Move construct from internal field, with specified patch type
        GeometricField
        (
            const IOobject& io,
            const Mesh& mesh,
            const dimensionSet& ds,
            Field<Type>&& iField,
            const word& patchFieldType = PatchField<Type>::calculatedType()
        );

        //- Copy construct from components
        GeometricField
        (
            const IOobject& io,
            const Mesh& mesh,
            const dimensionSet& ds,
            const Field<Type>& iField,
            const PtrList<PatchField<Type>>& ptfl
        );

        //- Construct and read given IOobject
        GeometricField
        (
            const IOobject& io,
            const Mesh& mesh,
            const bool readOldTime = true
        );

        //- Construct from dictionary
        GeometricField
        (
            const IOobject& io,
            const Mesh& mesh,
            const dictionary& dict
        );

        //- Copy construct
        GeometricField
        (
            const GeometricField<Type, PatchField, GeoMesh>& gf
        );

        //- Construct from tmp\<GeometricField\> deleting argument
        GeometricField
        (
            const tmp<GeometricField<Type, PatchField, GeoMesh>>& tgf
        );

        //- Construct as copy resetting IO parameters
        GeometricField
        (
            const IOobject& io,
            const GeometricField<Type, PatchField, GeoMesh>& gf
        );

        //- Construct as copy of tmp<GeometricField> resetting IO parameters
        GeometricField
        (
            const IOobject& io,
            const tmp<GeometricField<Type, PatchField, GeoMesh>>& tgf
        );

        //- Copy construct with a new name
        GeometricField
        (
            const word& newName,
            const GeometricField<Type, PatchField, GeoMesh>& gf
        );

        //- Construct with a new name from tmp\<GeometricField\>
        GeometricField
        (
            const word& newName,
            const tmp<GeometricField<Type, PatchField, GeoMesh>>& tgf
        );

        //- Construct as copy resetting IO parameters and patch type
        GeometricField
        (
            const IOobject& io,
            const GeometricField<Type, PatchField, GeoMesh>& gf,
            const word& patchFieldType
        );

        //- Construct as copy resetting IO parameters and boundary type
        //- for selected patchIDs
        GeometricField
        (
            const IOobject& io,
            const GeometricField<Type, PatchField, GeoMesh>& gf,
            const labelList& patchIDs,
            const word& patchFieldType
        );

        //- Construct as copy resetting IO parameters and boundary types
        GeometricField
        (
            const IOobject& io,
            const GeometricField<Type, PatchField, GeoMesh>& gf,
            const wordList& patchFieldTypes,
            const wordList& actualPatchTypes = wordList()
        );

        //- Construct as copy resetting IO parameters and boundary types
        GeometricField
        (
            const IOobject& io,
            const tmp<GeometricField<Type, PatchField, GeoMesh>>& tgf,
            const wordList& patchFieldTypes,
            const wordList& actualPatchTypes = wordList()
        );

        //- Clone
        tmp<GeometricField<Type, PatchField, GeoMesh>> clone() const;


    // Static Constructors

        //- Return tmp field from name, mesh, dimensions and patch type.
        //  The field is NO_READ, NO_WRITE, unregistered and uses the
        //  current timeName from the mesh registry
        static tmp<GeometricField<Type, PatchField, GeoMesh>> New
        (
            const word& name,
            const Mesh& mesh,
            const dimensionSet& ds,
            const word& patchFieldType = PatchField<Type>::calculatedType()
        );

        //- Return tmp field from name, mesh, dimensions,
        //- copy of internal field, with specified patch type.
        //  The field is NO_READ, NO_WRITE, unregistered and uses the
        //  current timeName from the mesh registry
        static tmp<GeometricField<Type, PatchField, GeoMesh>> New
        (
            const word& name,
            const Mesh& mesh,
            const dimensionSet& ds,
            const Field<Type>& iField,
            const word& patchFieldType = PatchField<Type>::calculatedType()
        );

        //- Return tmp field from name, mesh, dimensions,
        //- moved internal field contents, with specified patch type.
        //  The field is NO_READ, NO_WRITE, unregistered and uses the
        //  current timeName from the mesh registry
        static tmp<GeometricField<Type, PatchField, GeoMesh>> New
        (
            const word& name,
            const Mesh& mesh,
            const dimensionSet& ds,
            Field<Type>&& iField,
            const word& patchFieldType = PatchField<Type>::calculatedType()
        );

        //- Return tmp field from name, mesh, dimensioned\<Type\>
        //- and patch type.
        //  The field is NO_READ, NO_WRITE, unregistered and uses the
        //  current timeName from the mesh registry
        static tmp<GeometricField<Type, PatchField, GeoMesh>> New
        (
            const word& name,
            const Mesh& mesh,
            const dimensioned<Type>& dt,
            const word& patchFieldType = PatchField<Type>::calculatedType()
        );

        //- Return tmp field from name, mesh, dimensioned\<Type\>
        //- and patch types.
        //  The field is NO_READ, NO_WRITE, unregistered and uses the
        //  current timeName from the mesh registry
        static tmp<GeometricField<Type, PatchField, GeoMesh>> New
        (
            const word& name,
            const Mesh& mesh,
            const dimensioned<Type>& dt,
            const wordList& patchFieldTypes,
            const wordList& actualPatchTypes = wordList()
        );

        //- Return renamed tmp field
        //  The field is NO_READ, NO_WRITE, unregistered and uses the
        //  current timeName from the mesh registry
        static tmp<GeometricField<Type, PatchField, GeoMesh>> New
        (
            const word& newName,
            const tmp<GeometricField<Type, PatchField, GeoMesh>>& tgf
        );

        //- Rename tmp field and reset patch field types and return
        //  The field is NO_READ, NO_WRITE, unregistered and uses the
        //  current timeName from the mesh registry
        static tmp<GeometricField<Type, PatchField, GeoMesh>> New
        (
            const word& newName,
            const tmp<GeometricField<Type, PatchField, GeoMesh>>& tgf,
            const wordList& patchFieldTypes,
            const wordList& actualPatchTypes = wordList()
        );


    //- Destructor
    virtual ~GeometricField();


    // Member Functions

        //- Return a reference to the dimensioned internal field
        //  \param updateAccessTime update event counter and check
        //      old-time fields
        //
        //  \note Should avoid using updateAccessTime = true within loops.
        Internal& ref(const bool updateAccessTime = true);

        //- Return a const-reference to the dimensioned internal field
        inline const Internal& internalField() const;

        //- Return a const-reference to the dimensioned internal field
        //- of a "vol" field.
        //  Useful in the formulation of source-terms for FV equations
        //
        //  \note definition in finiteVolume/fields/volFields/volFieldsI.H
        inline const Internal& v() const;

        //- Return a reference to the internal field
        //  \param updateAccessTime update event counter and check
        //      old-time fields
        //
        //  \note Should avoid using updateAccessTime = true within loops.
        typename Internal::FieldType& primitiveFieldRef
        (
            const bool updateAccessTime = true
        );

        //- Return a const-reference to the  internal field
        inline const typename Internal::FieldType& primitiveField() const;

        //- Return a reference to the boundary field
        //  \param updateAccessTime update event counter and check
        //      old-time fields
        //
        //  \note Should avoid using updateAccessTime = true within loops.
        Boundary& boundaryFieldRef(const bool updateAccessTime = true);

        //- Return const-reference to the boundary field
        inline const Boundary& boundaryField() const;

        //- Return the time index of the field
        inline label timeIndex() const;

        //- Return the time index of the field
        inline label& timeIndex();

        //- Store the old-time fields
        void storeOldTimes() const;

        //- Store the old-time field
        void storeOldTime() const;

        //- Return the number of old time fields stored
        label nOldTimes() const;

        //- Return old time field
        const GeometricField<Type, PatchField, GeoMesh>& oldTime() const;

        //- Return non-const old time field
        //  (Not a good idea but it is used for sub-cycling)
        GeometricField<Type, PatchField, GeoMesh>& oldTime();

        //- Store the field as the previous iteration value
        void storePrevIter() const;

        //- Return previous iteration field
        const GeometricField<Type, PatchField, GeoMesh>& prevIter() const;

        //- Correct boundary field
        void correctBoundaryConditions();

        //- Does the field need a reference level for solution
        bool needReference() const;

        //- Return a component of the field
        tmp<GeometricField<cmptType, PatchField, GeoMesh>> component
        (
            const direction
        ) const;

        //- WriteData member function required by regIOobject
        bool writeData(Ostream&) const;

        //- Return transpose (only if it is a tensor field)
        tmp<GeometricField<Type, PatchField, GeoMesh>> T() const;

        //- Relax field (for steady-state solution).
        //  alpha = 1 : no relaxation
        //  alpha < 1 : relaxation
        //  alpha = 0 : do nothing
        void relax(const scalar alpha);

        //- Relax field (for steady-state solution).
        //  alpha is read from controlDict
        void relax();

        //- Select the final iteration parameters if `final' is true
        //  by returning the field name + "Final"
        //  otherwise the standard parameters by returning the field name
        word select(bool final) const;

        //- Helper function to write the min and max to an Ostream
        void writeMinMax(Ostream& os) const;


    // Member Function *this Operators

        //- Negate the field inplace
        void negate();

        //- Replace specified field component with content from another field
        void replace
        (
            const direction d,
            const GeometricField<cmptType, PatchField, GeoMesh>& gcf
        );

        //- Replace specified field component with specified value
        void replace
        (
            const direction d,
            const dimensioned<cmptType>& ds
        );

        //- Use the minimum of the field and specified value
        //  This sets the \em ceiling on the field values
        void min(const dimensioned<Type>& dt);

        //- Use the maximum of the field and specified value
        //  This sets the \em floor on the field values
        void max(const dimensioned<Type>& dt);

        //- Clip the field to be bounded within the specified range
        void clip(const dimensioned<MinMax<Type>>& range);

        //- Clip the field to be bounded within the specified range
        void clip
        (
            const dimensioned<Type>& minVal,
            const dimensioned<Type>& maxVal
        );

        //- Deprecated(2019-01) identical to clip()
        //  \deprecated(2019-01) identical to clip()
        FOAM_DEPRECATED_FOR(2019-01, "clip() method")
        void maxMin
        (
            const dimensioned<Type>& minVal,
            const dimensioned<Type>& maxVal
        );


    // Member Operators

        //- Return a const-reference to the dimensioned internal field
        //  Useful in the formulation of source-terms for FV equations
        inline const Internal& operator()() const;

        void operator=(const GeometricField<Type, PatchField, GeoMesh>&);
        void operator=(const tmp<GeometricField<Type, PatchField, GeoMesh>>&);
        void operator=(const dimensioned<Type>&);

        void operator==(const tmp<GeometricField<Type, PatchField, GeoMesh>>&);
        void operator==(const dimensioned<Type>&);

        void operator+=(const GeometricField<Type, PatchField, GeoMesh>&);
        void operator+=(const tmp<GeometricField<Type, PatchField, GeoMesh>>&);

        void operator-=(const GeometricField<Type, PatchField, GeoMesh>&);
        void operator-=(const tmp<GeometricField<Type, PatchField, GeoMesh>>&);

        void operator*=(const GeometricField<scalar, PatchField, GeoMesh>&);
        void operator*=(const tmp<GeometricField<scalar,PatchField,GeoMesh>>&);

        void operator/=(const GeometricField<scalar, PatchField, GeoMesh>&);
        void operator/=(const tmp<GeometricField<scalar,PatchField,GeoMesh>>&);

        void operator+=(const dimensioned<Type>&);
        void operator-=(const dimensioned<Type>&);

        void operator*=(const dimensioned<scalar>&);
        void operator/=(const dimensioned<scalar>&);


    // Ostream operators

        friend Ostream& operator<< <Type, PatchField, GeoMesh>
        (
            Ostream&,
            const GeometricField<Type, PatchField, GeoMesh>&
        );

        friend Ostream& operator<< <Type, PatchField, GeoMesh>
        (
            Ostream&,
            const tmp<GeometricField<Type, PatchField, GeoMesh>>&
        );
};


template<class Type, template<class> class PatchField, class GeoMesh>
Ostream& operator<<
(
    Ostream&,
    const typename GeometricField<Type, PatchField, GeoMesh>::
    Boundary&
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "GeometricFieldI.H"

#ifdef NoRepository
    #include "GeometricField.C"
#endif

#include "GeometricFieldFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
