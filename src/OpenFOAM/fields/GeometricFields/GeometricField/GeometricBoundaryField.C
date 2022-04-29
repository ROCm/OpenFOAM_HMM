/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017,2022 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "GeometricBoundaryField.H"
#include "globalMeshData.H"
#include "cyclicPolyPatch.H"
#include "emptyPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::readField
(
    const DimensionedField<Type, GeoMesh>& field,
    const dictionary& dict
)
{
    ///if (GeometricField<Type, PatchField, GeoMesh::debug)
    ///{
    ///    InfoInFunction << nl;
    ///}

    // Clear the boundary field if already initialised
    this->clear();

    this->resize(bmesh_.size());

    label nUnset = this->size();

    // 1. Handle explicit patch names. Note that there can be only one explicit
    //    patch name since is key of dictionary.

    for (const entry& dEntry : dict)
    {
        if (dEntry.isDict() && dEntry.keyword().isLiteral())
        {
            const label patchi = bmesh_.findPatchID(dEntry.keyword());

            if (patchi != -1)
            {
                this->set
                (
                    patchi,
                    PatchField<Type>::New
                    (
                        bmesh_[patchi],
                        field,
                        dEntry.dict()
                    )
                );
                nUnset--;
            }
        }
    }

    if (nUnset == 0)
    {
        return;
    }


    // 2. Patch-groups. (using non-wild card entries of dictionaries)
    // (patchnames already matched above)
    // Note: in reverse order of entries in the dictionary (last
    // patchGroups wins). This is so it is consistent with dictionary wildcard
    // behaviour
    for (auto iter = dict.crbegin(); iter != dict.crend(); ++iter)
    {
        const entry& dEntry = *iter;

        if (dEntry.isDict() && dEntry.keyword().isLiteral())
        {
            const labelList patchIds =
                bmesh_.indices(dEntry.keyword(), true); // use patchGroups

            for (const label patchi : patchIds)
            {
                if (!this->set(patchi))
                {
                    this->set
                    (
                        patchi,
                        PatchField<Type>::New
                        (
                            bmesh_[patchi],
                            field,
                            dEntry.dict()
                        )
                    );
                }
            }
        }
    }


    // 3. Wildcard patch overrides
    forAll(bmesh_, patchi)
    {
        if (!this->set(patchi))
        {
            if (bmesh_[patchi].type() == emptyPolyPatch::typeName)
            {
                this->set
                (
                    patchi,
                    PatchField<Type>::New
                    (
                        emptyPolyPatch::typeName,
                        bmesh_[patchi],
                        field
                    )
                );
            }
            else
            {
                bool found = dict.found(bmesh_[patchi].name());

                if (found)
                {
                    this->set
                    (
                        patchi,
                        PatchField<Type>::New
                        (
                            bmesh_[patchi],
                            field,
                            dict.subDict(bmesh_[patchi].name())
                        )
                    );
                }
            }
        }
    }


    // Check for any unset patches
    forAll(bmesh_, patchi)
    {
        if (!this->set(patchi))
        {
            if (bmesh_[patchi].type() == cyclicPolyPatch::typeName)
            {
                FatalIOErrorInFunction(dict)
                    << "Cannot find patchField entry for cyclic "
                    << bmesh_[patchi].name() << endl
                    << "Is your field uptodate with split cyclics?" << endl
                    << "Run foamUpgradeCyclics to convert mesh and fields"
                    << " to split cyclics." << exit(FatalIOError);
            }
            else
            {
                FatalIOErrorInFunction(dict)
                    << "Cannot find patchField entry for "
                    << bmesh_[patchi].name() << exit(FatalIOError);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::GeometricBoundaryField
(
    const BoundaryMesh& bmesh
)
:
    FieldField<PatchField, Type>(bmesh.size()),
    bmesh_(bmesh)
{}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::GeometricBoundaryField
(
    const BoundaryMesh& bmesh,
    const DimensionedField<Type, GeoMesh>& field,
    const word& patchFieldType
)
:
    FieldField<PatchField, Type>(bmesh.size()),
    bmesh_(bmesh)
{
    ///if (GeometricField<Type, PatchField, GeoMesh::debug)
    ///{
    ///    InfoInFunction << nl;
    ///}

    forAll(bmesh_, patchi)
    {
        this->set
        (
            patchi,
            PatchField<Type>::New
            (
                patchFieldType,
                bmesh_[patchi],
                field
            )
        );
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::GeometricBoundaryField
(
    const BoundaryMesh& bmesh,
    const DimensionedField<Type, GeoMesh>& field,
    const wordList& patchFieldTypes,
    const wordList& constraintTypes
)
:
    FieldField<PatchField, Type>(bmesh.size()),
    bmesh_(bmesh)
{
    ///if (GeometricField<Type, PatchField, GeoMesh::debug)
    ///{
    ///    InfoInFunction << nl;
    ///}

    if
    (
        patchFieldTypes.size() != this->size()
     || (constraintTypes.size() && (constraintTypes.size() != this->size()))
    )
    {
        FatalErrorInFunction
            << "Incorrect number of patch type specifications given" << nl
            << "    Number of patches in mesh = " << bmesh.size()
            << " number of patch type specifications = "
            << patchFieldTypes.size()
            << abort(FatalError);
    }

    if (constraintTypes.size())
    {
        forAll(bmesh_, patchi)
        {
            this->set
            (
                patchi,
                PatchField<Type>::New
                (
                    patchFieldTypes[patchi],
                    constraintTypes[patchi],
                    bmesh_[patchi],
                    field
                )
            );
        }
    }
    else
    {
        forAll(bmesh_, patchi)
        {
            this->set
            (
                patchi,
                PatchField<Type>::New
                (
                    patchFieldTypes[patchi],
                    bmesh_[patchi],
                    field
                )
            );
        }
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::GeometricBoundaryField
(
    const BoundaryMesh& bmesh,
    const DimensionedField<Type, GeoMesh>& field,
    const PtrList<PatchField<Type>>& ptfl
)
:
    FieldField<PatchField, Type>(bmesh.size()),
    bmesh_(bmesh)
{
    ///if (GeometricField<Type, PatchField, GeoMesh::debug)
    ///{
    ///    InfoInFunction << nl;
    ///}

    forAll(bmesh_, patchi)
    {
        this->set(patchi, ptfl[patchi].clone(field));
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::GeometricBoundaryField
(
    const DimensionedField<Type, GeoMesh>& field,
    const GeometricBoundaryField<Type, PatchField, GeoMesh>& btf
)
:
    FieldField<PatchField, Type>(btf.size()),
    bmesh_(btf.bmesh_)
{
    ///if (GeometricField<Type, PatchField, GeoMesh::debug)
    ///{
    ///    InfoInFunction << nl;
    ///}

    forAll(bmesh_, patchi)
    {
        this->set(patchi, btf[patchi].clone(field));
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::GeometricBoundaryField
(
    const DimensionedField<Type, GeoMesh>& field,
    const GeometricBoundaryField<Type, PatchField, GeoMesh>& btf,
    const labelList& patchIDs,
    const word& patchFieldType
)
:
    FieldField<PatchField, Type>(btf.size()),
    bmesh_(btf.bmesh_)
{
    ///if (GeometricField<Type, PatchField, GeoMesh::debug)
    ///{
    ///    InfoInFunction << nl;
    ///}

    for (const label patchi : patchIDs)
    {
        this->set
        (
            patchi,
            PatchField<Type>::New
            (
                patchFieldType,
                bmesh_[patchi],
                field
            )
        );
    }

    forAll(bmesh_, patchi)
    {
        if (!this->set(patchi))
        {
            this->set(patchi, btf[patchi].clone(field));
        }
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::GeometricBoundaryField
(
    const GeometricBoundaryField<Type, PatchField, GeoMesh>& btf
)
:
    FieldField<PatchField, Type>(btf),
    bmesh_(btf.bmesh_)
{
    ///if (GeometricField<Type, PatchField, GeoMesh::debug)
    ///{
    ///    InfoInFunction << nl;
    ///}
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::GeometricBoundaryField
(
    const BoundaryMesh& bmesh,
    const DimensionedField<Type, GeoMesh>& field,
    const dictionary& dict
)
:
    FieldField<PatchField, Type>(bmesh.size()),
    bmesh_(bmesh)
{
    readField(field, dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::updateCoeffs()
{
    ///if (GeometricField<Type, PatchField, GeoMesh::debug)
    ///{
    ///    InfoInFunction << nl;
    ///}

    for (auto& pfld : *this)
    {
        pfld.updateCoeffs();
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::evaluate()
{
    ///if (GeometricField<Type, PatchField, GeoMesh::debug)
    ///{
    ///    InfoInFunction << nl;
    ///}

    const UPstream::commsTypes commsType(UPstream::defaultCommsType);

    if
    (
        commsType == UPstream::commsTypes::blocking
     || commsType == UPstream::commsTypes::nonBlocking
    )
    {
        const label startOfRequests = UPstream::nRequests();

        for (auto& pfld : *this)
        {
            pfld.initEvaluate(commsType);
        }

        // Wait for outstanding requests
        if
        (
            UPstream::parRun()
         && commsType == Pstream::commsTypes::nonBlocking
        )
        {
            UPstream::waitRequests(startOfRequests);
        }

        for (auto& pfld : *this)
        {
            pfld.evaluate(commsType);
        }
    }
    else if (commsType == UPstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule =
            bmesh_.mesh().globalData().patchSchedule();

        for (const auto& schedEval : patchSchedule)
        {
            const label patchi = schedEval.patch;
            auto& pfld = (*this)[patchi];

            if (schedEval.init)
            {
                pfld.initEvaluate(commsType);
            }
            else
            {
                pfld.evaluate(commsType);
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported communications type "
            << UPstream::commsTypeNames[commsType]
            << exit(FatalError);
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
template<class CoupledPatchType>
void Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::evaluateCoupled()
{
    ///if (GeometricField<Type, PatchField, GeoMesh::debug)
    ///{
    ///    InfoInFunction << nl;
    ///}

    const UPstream::commsTypes commsType(UPstream::defaultCommsType);

    if
    (
        commsType == UPstream::commsTypes::blocking
     || commsType == UPstream::commsTypes::nonBlocking
    )
    {
        const label startOfRequests = UPstream::nRequests();

        for (auto& pfld : *this)
        {
            const auto* cpp = isA<CoupledPatchType>(pfld.patch());

            if (cpp && cpp->coupled())
            {
                pfld.initEvaluate(commsType);
            }
        }

        // Wait for outstanding requests
        if
        (
            UPstream::parRun()
         && commsType == UPstream::commsTypes::nonBlocking
        )
        {
            UPstream::waitRequests(startOfRequests);
        }

        for (auto& pfld : *this)
        {
            const auto* cpp = isA<CoupledPatchType>(pfld.patch());

            if (cpp && cpp->coupled())
            {
                pfld.evaluate(commsType);
            }
        }
    }
    else if (commsType == UPstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule =
            bmesh_.mesh().globalData().patchSchedule();

        for (const auto& schedEval : patchSchedule)
        {
            const label patchi = schedEval.patch;
            auto& pfld = (*this)[patchi];

            const auto* cpp = isA<CoupledPatchType>(pfld.patch());

            if (cpp && cpp->coupled())
            {
                if (schedEval.init)
                {
                    pfld.initEvaluate(commsType);
                }
                else
                {
                    pfld.evaluate(commsType);
                }
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported communications type "
            << UPstream::commsTypeNames[commsType]
            << exit(FatalError);
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::wordList
Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::types() const
{
    const FieldField<PatchField, Type>& pff = *this;

    wordList list(pff.size());

    forAll(pff, patchi)
    {
        list[patchi] = pff[patchi].type();
    }

    return list;
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>
Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::
boundaryInternalField() const
{
    GeometricBoundaryField<Type, PatchField, GeoMesh> result(*this);

    forAll(result, patchi)
    {
        result[patchi] == this->operator[](patchi).patchInternalField();
    }

    return result;
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::LduInterfaceFieldPtrsList<Type>
Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::interfaces() const
{
    LduInterfaceFieldPtrsList<Type> list(this->size());

    forAll(list, patchi)
    {
        const auto* lduPtr =
            isA<LduInterfaceField<Type>>(this->operator[](patchi));

        if (lduPtr)
        {
            list.set(patchi, lduPtr);
        }
    }

    return list;
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::lduInterfaceFieldPtrsList
Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::
scalarInterfaces() const
{
    lduInterfaceFieldPtrsList list(this->size());

    forAll(list, patchi)
    {
        const auto* lduPtr =
            isA<lduInterfaceField>(this->operator[](patchi));

        if (lduPtr)
        {
            list.set(patchi, lduPtr);
        }
    }

    return list;
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::writeEntry
(
    const word& keyword,
    Ostream& os
) const
{
    os.beginBlock(keyword);
    this->writeEntries(os);
    os.endBlock();

    os.check(FUNCTION_NAME);
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::writeEntries
(
    Ostream& os
) const
{
    for (const auto& pfld : *this)
    {
        os.beginBlock(pfld.patch().name());
        os << pfld;
        os.endBlock();
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::operator=
(
    const GeometricBoundaryField<Type, PatchField, GeoMesh>& bf
)
{
    FieldField<PatchField, Type>::operator=(bf);
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::operator=
(
    const FieldField<PatchField, Type>& bf
)
{
    FieldField<PatchField, Type>::operator=(bf);
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::operator=
(
    const Type& val
)
{
    FieldField<PatchField, Type>::operator=(val);
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::operator==
(
    const GeometricBoundaryField<Type, PatchField, GeoMesh>& bf
)
{
    forAll(*this, patchi)
    {
        this->operator[](patchi) == bf[patchi];
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::operator==
(
    const FieldField<PatchField, Type>& bf
)
{
    forAll(*this, patchi)
    {
        this->operator[](patchi) == bf[patchi];
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricBoundaryField<Type, PatchField, GeoMesh>::operator==
(
    const Type& val
)
{
    forAll(*this, patchi)
    {
        this->operator[](patchi) == val;
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const GeometricBoundaryField<Type, PatchField, GeoMesh>& bf
)
{
    os << static_cast<const FieldField<PatchField, Type>&>(bf);
    return os;
}


// ************************************************************************* //
