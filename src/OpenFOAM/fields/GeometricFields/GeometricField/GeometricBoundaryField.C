/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "emptyPolyPatch.H"
#include "commSchedule.H"
#include "globalMeshData.H"
#include "cyclicPolyPatch.H"

template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::Boundary::readField
(
    const DimensionedField<Type, GeoMesh>& field,
    const dictionary& dict
)
{
    DebugInFunction << nl;

    // Clear the boundary field if already initialised
    this->clear();

    this->setSize(bmesh_.size());

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
Foam::GeometricField<Type, PatchField, GeoMesh>::Boundary::Boundary
(
    const BoundaryMesh& bmesh
)
:
    FieldField<PatchField, Type>(bmesh.size()),
    bmesh_(bmesh)
{}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::Boundary::Boundary
(
    const BoundaryMesh& bmesh,
    const DimensionedField<Type, GeoMesh>& field,
    const word& patchFieldType
)
:
    FieldField<PatchField, Type>(bmesh.size()),
    bmesh_(bmesh)
{
    DebugInFunction << nl;

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
Foam::GeometricField<Type, PatchField, GeoMesh>::Boundary::Boundary
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
    DebugInFunction << nl;

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
Foam::GeometricField<Type, PatchField, GeoMesh>::Boundary::Boundary
(
    const BoundaryMesh& bmesh,
    const DimensionedField<Type, GeoMesh>& field,
    const PtrList<PatchField<Type>>& ptfl
)
:
    FieldField<PatchField, Type>(bmesh.size()),
    bmesh_(bmesh)
{
    DebugInFunction << nl;

    forAll(bmesh_, patchi)
    {
        this->set(patchi, ptfl[patchi].clone(field));
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::Boundary::Boundary
(
    const DimensionedField<Type, GeoMesh>& field,
    const typename GeometricField<Type, PatchField, GeoMesh>::Boundary& btf
)
:
    FieldField<PatchField, Type>(btf.size()),
    bmesh_(btf.bmesh_)
{
    DebugInFunction << nl;

    forAll(bmesh_, patchi)
    {
        this->set(patchi, btf[patchi].clone(field));
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::Boundary::Boundary
(
    const DimensionedField<Type, GeoMesh>& field,
    const typename GeometricField<Type, PatchField, GeoMesh>::Boundary& btf,
    const labelList& patchIDs,
    const word& patchFieldType
)
:
    FieldField<PatchField, Type>(btf.size()),
    bmesh_(btf.bmesh_)
{
    DebugInFunction << nl;

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
Foam::GeometricField<Type, PatchField, GeoMesh>::Boundary::Boundary
(
    const typename GeometricField<Type, PatchField, GeoMesh>::Boundary& btf
)
:
    FieldField<PatchField, Type>(btf),
    bmesh_(btf.bmesh_)
{
    DebugInFunction << nl;
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::GeometricField<Type, PatchField, GeoMesh>::Boundary::Boundary
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
void Foam::GeometricField<Type, PatchField, GeoMesh>::Boundary::updateCoeffs()
{
    DebugInFunction << nl;

    forAll(*this, patchi)
    {
        this->operator[](patchi).updateCoeffs();
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::Boundary::evaluate()
{
    DebugInFunction << nl;

    if
    (
        Pstream::defaultCommsType == Pstream::commsTypes::blocking
     || Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
    )
    {
        const label nReq = Pstream::nRequests();

        forAll(*this, patchi)
        {
            this->operator[](patchi).initEvaluate(Pstream::defaultCommsType);
        }

        // Block for any outstanding requests
        if
        (
            Pstream::parRun()
         && Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
        )
        {
            Pstream::waitRequests(nReq);
        }

        forAll(*this, patchi)
        {
            this->operator[](patchi).evaluate(Pstream::defaultCommsType);
        }
    }
    else if (Pstream::defaultCommsType == Pstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule =
            bmesh_.mesh().globalData().patchSchedule();

        forAll(patchSchedule, patchEvali)
        {
            if (patchSchedule[patchEvali].init)
            {
                this->operator[](patchSchedule[patchEvali].patch)
                    .initEvaluate(Pstream::commsTypes::scheduled);
            }
            else
            {
                this->operator[](patchSchedule[patchEvali].patch)
                    .evaluate(Pstream::commsTypes::scheduled);
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unsupported communications type "
            << Pstream::commsTypeNames[Pstream::defaultCommsType]
            << exit(FatalError);
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::wordList
Foam::GeometricField<Type, PatchField, GeoMesh>::Boundary::types() const
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
typename Foam::GeometricField<Type, PatchField, GeoMesh>::Boundary
Foam::GeometricField<Type, PatchField, GeoMesh>::Boundary::
boundaryInternalField() const
{
    typename GeometricField<Type, PatchField, GeoMesh>::Boundary
        BoundaryInternalField(*this);

    forAll(BoundaryInternalField, patchi)
    {
        BoundaryInternalField[patchi] ==
            this->operator[](patchi).patchInternalField();
    }

    return BoundaryInternalField;
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::LduInterfaceFieldPtrsList<Type>
Foam::GeometricField<Type, PatchField, GeoMesh>::Boundary::interfaces() const
{
    LduInterfaceFieldPtrsList<Type> list(this->size());

    forAll(list, patchi)
    {
        const LduInterfaceField<Type>* lduPtr =
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
Foam::GeometricField<Type, PatchField, GeoMesh>::Boundary::
scalarInterfaces() const
{
    lduInterfaceFieldPtrsList list(this->size());

    forAll(list, patchi)
    {
        const lduInterfaceField* lduPtr =
            isA<lduInterfaceField>(this->operator[](patchi));

        if (lduPtr)
        {
            list.set(patchi, lduPtr);
        }
    }

    return list;
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::Boundary::writeEntry
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
void Foam::GeometricField<Type, PatchField, GeoMesh>::Boundary::writeEntries
(
    Ostream& os
) const
{
    forAll(*this, patchi)
    {
        os.beginBlock(this->operator[](patchi).patch().name());
        os  << this->operator[](patchi);
        os.endBlock();
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::Boundary::operator=
(
    const typename GeometricField<Type, PatchField, GeoMesh>::
    Boundary& bf
)
{
    FieldField<PatchField, Type>::operator=(bf);
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::Boundary::operator=
(
    const FieldField<PatchField, Type>& ptff
)
{
    FieldField<PatchField, Type>::operator=(ptff);
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::Boundary::operator=
(
    const Type& t
)
{
    FieldField<PatchField, Type>::operator=(t);
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::Boundary::operator==
(
    const typename GeometricField<Type, PatchField, GeoMesh>::
    Boundary& bf
)
{
    forAll(*this, patchi)
    {
        this->operator[](patchi) == bf[patchi];
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::Boundary::operator==
(
    const FieldField<PatchField, Type>& ptff
)
{
    forAll(*this, patchi)
    {
        this->operator[](patchi) == ptff[patchi];
    }
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::GeometricField<Type, PatchField, GeoMesh>::Boundary::operator==
(
    const Type& t
)
{
    forAll(*this, patchi)
    {
        this->operator[](patchi) == t;
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const typename GeometricField<Type, PatchField, GeoMesh>::
    Boundary& bf
)
{
    os << static_cast<const FieldField<PatchField, Type>&>(bf);
    return os;
}


// ************************************************************************* //
