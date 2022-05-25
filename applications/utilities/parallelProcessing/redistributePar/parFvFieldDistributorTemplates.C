/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
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

#include "parFvFieldDistributor.H"
#include "Time.H"
#include "PtrList.H"
#include "fvPatchFields.H"
#include "emptyFvPatch.H"
#include "emptyFvPatchField.H"
#include "emptyFvsPatchField.H"
#include "IOobjectList.H"
#include "mapDistributePolyMesh.H"
#include "processorFvPatch.H"

#include "distributedFieldMapper.H"
#include "distributedFvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::DimensionedField<Type, Foam::volMesh>>
Foam::parFvFieldDistributor::distributeField
(
    const DimensionedField<Type, volMesh>& fld
) const
{
    // Create internalField by remote mapping

    distributedFieldMapper mapper
    (
        labelUList::null(),
        distMap_.cellMap()
    );

    auto tfield = tmp<DimensionedField<Type, volMesh>>::New
    (
        IOobject
        (
            fld.name(),
            tgtMesh_.time().timeName(),
            fld.local(),
            tgtMesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        tgtMesh_,
        fld.dimensions(),
        Field<Type>(fld, mapper)
    );

    tfield.ref().oriented() = fld.oriented();

    return tfield;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::parFvFieldDistributor::distributeField
(
    const GeometricField<Type, fvPatchField, volMesh>& fld
) const
{
    // Create internalField by remote mapping
    distributedFieldMapper mapper
    (
        labelUList::null(),
        distMap_.cellMap()
    );

    DimensionedField<Type, volMesh> internalField
    (
        IOobject
        (
            fld.name(),
            tgtMesh_.time().timeName(),
            fld.local(),
            tgtMesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        tgtMesh_,
        fld.dimensions(),
        Field<Type>(fld.internalField(), mapper)
    );

    internalField.oriented() = fld.oriented();


    // Create patchFields by remote mapping
    // Note: patchFields still on source mesh, not target mesh

    PtrList<fvPatchField<Type>> oldPatchFields(fld.mesh().boundary().size());

    const auto& bfld = fld.boundaryField();

    forAll(bfld, patchi)
    {
        if (patchFaceMaps_.set(patchi))
        {
            // Clone local patch field
            oldPatchFields.set(patchi, bfld[patchi].clone());

            distributedFvPatchFieldMapper mapper
            (
                labelUList::null(),
                patchFaceMaps_[patchi]
            );

            // Map into local copy
            oldPatchFields[patchi].autoMap(mapper);
        }
    }


    // Clone the oldPatchFields onto the target patches. This is just to reset
    // the reference to the patch, size and content stay the same.

    PtrList<fvPatchField<Type>> newPatchFields(tgtMesh_.boundary().size());

    forAll(oldPatchFields, patchi)
    {
        if (oldPatchFields.set(patchi))
        {
            const auto& pfld = oldPatchFields[patchi];

            labelList dummyMap(identity(pfld.size()));
            directFvPatchFieldMapper dummyMapper(dummyMap);

            newPatchFields.set
            (
                patchi,
                fvPatchField<Type>::New
                (
                    pfld,
                    tgtMesh_.boundary()[patchi],
                    DimensionedField<Type, volMesh>::null(),
                    dummyMapper
                )
            );
        }
    }

    // Add some empty patches on remaining patches
    // (... probably processor patches)

    forAll(newPatchFields, patchi)
    {
        if (!newPatchFields.set(patchi))
        {
            newPatchFields.set
            (
                patchi,
                fvPatchField<Type>::New
                (
                    emptyFvPatchField<Type>::typeName,
                    tgtMesh_.boundary()[patchi],
                    DimensionedField<Type, volMesh>::null()
                )
            );
        }
    }

    // Return geometric field

    return tmp<GeometricField<Type, fvPatchField, volMesh>>::New
    (
        std::move(internalField),
        newPatchFields
    );
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::parFvFieldDistributor::distributeField
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& fld
) const
{
    // Create internalField by remote mapping
    distributedFieldMapper mapper
    (
        labelUList::null(),
        distMap_.faceMap()
    );


    Field<Type> primitiveField;
    {
        // Create flat field of internalField + all patch fields
        Field<Type> flatFld(fld.mesh().nFaces(), Type(Zero));
        SubList<Type>(flatFld, fld.internalField().size())
            = fld.internalField();

        for (const fvsPatchField<Type>& fvp : fld.boundaryField())
        {
            SubList<Type>(flatFld, fvp.size(), fvp.patch().start()) = fvp;
        }

        // Map all faces
        primitiveField = Field<Type>(flatFld, mapper, fld.oriented()());

        // Trim to internal faces (note: could also have special mapper)
        primitiveField.resize
        (
            min
            (
                primitiveField.size(),
                tgtMesh_.nInternalFaces()
            )
        );
    }


    DimensionedField<Type, surfaceMesh> internalField
    (
        IOobject
        (
            fld.name(),
            tgtMesh_.time().timeName(),
            fld.local(),
            tgtMesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        tgtMesh_,
        fld.dimensions(),
        std::move(primitiveField)
    );

    internalField.oriented() = fld.oriented();


    // Create patchFields by remote mapping
    // Note: patchFields still on source mesh, not target mesh

    PtrList<fvsPatchField<Type>> oldPatchFields(fld.mesh().boundary().size());

    const auto& bfld = fld.boundaryField();

    forAll(bfld, patchi)
    {
        if (patchFaceMaps_.set(patchi))
        {
            // Clone local patch field
            oldPatchFields.set(patchi, bfld[patchi].clone());

            distributedFvPatchFieldMapper mapper
            (
                labelUList::null(),
                patchFaceMaps_[patchi]
            );

            // Map into local copy
            oldPatchFields[patchi].autoMap(mapper);
        }
    }


    PtrList<fvsPatchField<Type>> newPatchFields(tgtMesh_.boundary().size());

    // Clone the patchFields onto the base patches. This is just to reset
    // the reference to the patch, size and content stay the same.
    forAll(oldPatchFields, patchi)
    {
        if (oldPatchFields.set(patchi))
        {
            const fvsPatchField<Type>& pfld = oldPatchFields[patchi];

            labelList dummyMap(identity(pfld.size()));
            directFvPatchFieldMapper dummyMapper(dummyMap);

            newPatchFields.set
            (
                patchi,
                fvsPatchField<Type>::New
                (
                    pfld,
                    tgtMesh_.boundary()[patchi],
                    DimensionedField<Type, surfaceMesh>::null(),
                    dummyMapper
                )
            );
        }
    }

    // Add some empty patches on remaining patches
    // (... probably processor patches)
    forAll(newPatchFields, patchi)
    {
        if (!newPatchFields.set(patchi))
        {
            newPatchFields.set
            (
                patchi,
                fvsPatchField<Type>::New
                (
                    emptyFvsPatchField<Type>::typeName,
                    tgtMesh_.boundary()[patchi],
                    DimensionedField<Type, surfaceMesh>::null()
                )
            );
        }
    }


    // Return geometric field
    return tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>::New
    (
        std::move(internalField),
        newPatchFields
    );
}


template<class Type>
Foam::tmp<Foam::DimensionedField<Type, Foam::volMesh>>
Foam::parFvFieldDistributor::distributeInternalField
(
    const IOobject& fieldObject
) const
{
    // Read field
    DimensionedField<Type, volMesh> fld
    (
        fieldObject,
        srcMesh_
    );

    // Distribute
    return distributeField(fld);
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::parFvFieldDistributor::distributeVolumeField
(
    const IOobject& fieldObject
) const
{
    // Read field
    GeometricField<Type, fvPatchField, volMesh> fld
    (
        fieldObject,
        srcMesh_
    );

    // Distribute
    return distributeField(fld);
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::parFvFieldDistributor::distributeSurfaceField
(
    const IOobject& fieldObject
) const
{
    // Read field
    GeometricField<Type, fvsPatchField, surfaceMesh> fld
    (
        fieldObject,
        srcMesh_
    );

    // Distribute
    return distributeField(fld);
}


template<class Type>
Foam::label Foam::parFvFieldDistributor::distributeInternalFields
(
    const IOobjectList& objects,
    const wordRes& selectedFields
) const
{
    typedef DimensionedField<Type, volMesh> fieldType;

    // Available fields, sorted order
    const wordList fieldNames =
    (
        selectedFields.empty()
      ? objects.sortedNames<fieldType>()
      : objects.sortedNames<fieldType>(selectedFields)
    );

    label nFields = 0;
    for (const word& fieldName : fieldNames)
    {
        if ("cellDist" == fieldName)
        {
            // There is an odd chance this is an internal field
            continue;
        }
        if (verbose_)
        {
            if (!nFields)
            {
                Info<< "    Reconstructing "
                    << fieldType::typeName << "s\n" << nl;
            }
            Info<< "        " << fieldName << nl;
        }
        ++nFields;

        tmp<fieldType> tfld
        (
            distributeInternalField<Type>(*(objects[fieldName]))
        );
        if (isWriteProc_)
        {
            tfld().write();
        }
    }

    if (nFields && verbose_) Info<< endl;
    return nFields;
}


template<class Type>
Foam::label Foam::parFvFieldDistributor::distributeVolumeFields
(
    const IOobjectList& objects,
    const wordRes& selectedFields
) const
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    // Available fields, sorted order
    const wordList fieldNames =
    (
        selectedFields.empty()
      ? objects.sortedNames<fieldType>()
      : objects.sortedNames<fieldType>(selectedFields)
    );

    label nFields = 0;
    for (const word& fieldName : fieldNames)
    {
        if ("cellDist" == fieldName)
        {
            continue;
        }
        if (verbose_)
        {
            if (!nFields)
            {
                Info<< "    Reconstructing "
                    << fieldType::typeName << "s\n" << nl;
            }
            Info<< "        " << fieldName << nl;
        }
        ++nFields;

        tmp<fieldType> tfld
        (
            distributeVolumeField<Type>(*(objects[fieldName]))
        );
        if (isWriteProc_)
        {
            tfld().write();
        }
    }

    if (nFields && verbose_) Info<< endl;
    return nFields;
}


template<class Type>
Foam::label Foam::parFvFieldDistributor::distributeSurfaceFields
(
    const IOobjectList& objects,
    const wordRes& selectedFields
) const
{
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> fieldType;

    // Available fields, sorted order
    const wordList fieldNames =
    (
        selectedFields.empty()
      ? objects.sortedNames<fieldType>()
      : objects.sortedNames<fieldType>(selectedFields)
    );

    label nFields = 0;
    for (const word& fieldName : fieldNames)
    {
        if (verbose_)
        {
            if (!nFields)
            {
                Info<< "    Reconstructing "
                    << fieldType::typeName << "s\n" << nl;
            }
            Info<< "        " << fieldName << nl;
        }
        ++nFields;

        tmp<fieldType> tfld
        (
            distributeSurfaceField<Type>(*(objects[fieldName]))
        );
        if (isWriteProc_)
        {
            tfld().write();
        }
    }

    if (nFields && verbose_) Info<< endl;
    return nFields;
}


// ************************************************************************* //
