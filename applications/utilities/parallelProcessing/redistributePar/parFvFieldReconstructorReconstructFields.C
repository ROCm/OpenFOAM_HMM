/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
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

#include "parFvFieldReconstructor.H"
#include "Time.H"
#include "PtrList.H"
#include "fvPatchFields.H"
#include "emptyFvPatch.H"
#include "emptyFvPatchField.H"
#include "emptyFvsPatchField.H"
#include "IOobjectList.H"
#include "mapDistributePolyMesh.H"
#include "processorFvPatch.H"

#include "directFvPatchFieldMapper.H"
#include "distributedUnallocatedDirectFieldMapper.H"
#include "distributedUnallocatedDirectFvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::DimensionedField<Type, Foam::volMesh>>
Foam::parFvFieldReconstructor::reconstructFvVolumeInternalField
(
    const DimensionedField<Type, volMesh>& fld
) const
{
    distributedUnallocatedDirectFieldMapper mapper
    (
        labelUList::null(),
        distMap_.cellMap()
    );

    Field<Type> internalField(fld, mapper);

    // Construct a volField
    IOobject baseIO
    (
        fld.name(),
        baseMesh_.time().timeName(),
        fld.local(),
        baseMesh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    tmp<DimensionedField<Type, volMesh>> tfield
    (
        new DimensionedField<Type, volMesh>
        (
            baseIO,
            baseMesh_,
            fld.dimensions(),
            internalField
        )
    );

    tfield.ref().oriented() = fld.oriented();

    return tfield;
}


template<class Type>
Foam::tmp<Foam::DimensionedField<Type, Foam::volMesh>>
Foam::parFvFieldReconstructor::reconstructFvVolumeInternalField
(
    const IOobject& fieldIoObject
) const
{
    // Read the field
    DimensionedField<Type, volMesh> fld
    (
        fieldIoObject,
        procMesh_
    );

    // Distribute onto baseMesh
    return reconstructFvVolumeInternalField(fld);
}


// Reconstruct a field onto the baseMesh
template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::parFvFieldReconstructor::reconstructFvVolumeField
(
    const GeometricField<Type, fvPatchField, volMesh>& fld
) const
{
    // Create the internalField by remote mapping
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    distributedUnallocatedDirectFieldMapper mapper
    (
        labelUList::null(),
        distMap_.cellMap()
    );

    Field<Type> internalField(fld.internalField(), mapper);



    // Create the patchFields by remote mapping
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Note: patchFields still on mesh, not baseMesh

    PtrList<fvPatchField<Type>> patchFields(fld.mesh().boundary().size());

    const typename GeometricField<Type, fvPatchField, volMesh>::Boundary&
        bfld = fld.boundaryField();

    forAll(bfld, patchI)
    {
        if (patchFaceMaps_.set(patchI))
        {
            // Clone local patch field
            patchFields.set(patchI, bfld[patchI].clone());

            distributedUnallocatedDirectFvPatchFieldMapper mapper
            (
                labelUList::null(),
                patchFaceMaps_[patchI]
            );

            // Map into local copy
            patchFields[patchI].autoMap(mapper);
        }
    }


    PtrList<fvPatchField<Type>> basePatchFields
    (
        baseMesh_.boundary().size()
    );

    // Clone the patchFields onto the base patches. This is just to reset
    // the reference to the patch, size and content stay the same.
    forAll(patchFields, patchI)
    {
        if (patchFields.set(patchI))
        {
            const fvPatch& basePatch = baseMesh_.boundary()[patchI];

            const fvPatchField<Type>& pfld = patchFields[patchI];

            labelList dummyMap(identity(pfld.size()));
            directFvPatchFieldMapper dummyMapper(dummyMap);

            basePatchFields.set
            (
                patchI,
                fvPatchField<Type>::New
                (
                    pfld,
                    basePatch,
                    DimensionedField<Type, volMesh>::null(),
                    dummyMapper
                )
            );
        }
    }

    // Add some empty patches on remaining patches (tbd.probably processor
    // patches)
    forAll(basePatchFields, patchI)
    {
        if (patchI >= patchFields.size() || !patchFields.set(patchI))
        {
            basePatchFields.set
            (
                patchI,
                fvPatchField<Type>::New
                (
                    emptyFvPatchField<Type>::typeName,
                    baseMesh_.boundary()[patchI],
                    DimensionedField<Type, volMesh>::null()
                )
            );
        }
    }

    // Construct a volField
    IOobject baseIO
    (
        fld.name(),
        baseMesh_.time().timeName(),
        fld.local(),
        baseMesh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    tmp<GeometricField<Type, fvPatchField, volMesh>> tfield
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            baseIO,
            baseMesh_,
            fld.dimensions(),
            internalField,
            basePatchFields
        )
    );

    tfield.ref().oriented()= fld.oriented();

    return tfield;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh>>
Foam::parFvFieldReconstructor::reconstructFvVolumeField
(
    const IOobject& fieldIoObject
) const
{
    // Read the field
    GeometricField<Type, fvPatchField, volMesh> fld
    (
        fieldIoObject,
        procMesh_
    );

    // Distribute onto baseMesh
    return reconstructFvVolumeField(fld);
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::parFvFieldReconstructor::reconstructFvSurfaceField
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& fld
) const
{
    // Create the internalField by remote mapping
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    distributedUnallocatedDirectFieldMapper mapper
    (
        labelUList::null(),
        distMap_.faceMap()
    );

    // Create flat field of internalField + all patch fields
    Field<Type> flatFld(fld.mesh().nFaces(), Type(Zero));
    SubList<Type>(flatFld, fld.internalField().size()) = fld.internalField();
    forAll(fld.boundaryField(), patchI)
    {
        const fvsPatchField<Type>& fvp = fld.boundaryField()[patchI];

        SubList<Type>(flatFld, fvp.size(), fvp.patch().start()) = fvp;
    }

    // Map all faces
    Field<Type> internalField(flatFld, mapper, fld.oriented()());

    // Trim to internal faces (note: could also have special mapper)
    internalField.setSize
    (
        min
        (
            internalField.size(),
            baseMesh_.nInternalFaces()
        )
    );


    // Create the patchFields by remote mapping
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Note: patchFields still on mesh, not baseMesh

    PtrList<fvsPatchField<Type>> patchFields(fld.mesh().boundary().size());

    const typename GeometricField<Type, fvsPatchField, surfaceMesh>::Boundary&
        bfld = fld.boundaryField();

    forAll(bfld, patchI)
    {
        if (patchFaceMaps_.set(patchI))
        {
            // Clone local patch field
            patchFields.set(patchI, bfld[patchI].clone());

            distributedUnallocatedDirectFvPatchFieldMapper mapper
            (
                labelUList::null(),
                patchFaceMaps_[patchI]
            );

            // Map into local copy
            patchFields[patchI].autoMap(mapper);
        }
    }


    PtrList<fvsPatchField<Type>> basePatchFields
    (
        baseMesh_.boundary().size()
    );

    // Clone the patchFields onto the base patches. This is just to reset
    // the reference to the patch, size and content stay the same.
    forAll(patchFields, patchI)
    {
        if (patchFields.set(patchI))
        {
            const fvPatch& basePatch = baseMesh_.boundary()[patchI];

            const fvsPatchField<Type>& pfld = patchFields[patchI];

            labelList dummyMap(identity(pfld.size()));
            directFvPatchFieldMapper dummyMapper(dummyMap);

            basePatchFields.set
            (
                patchI,
                fvsPatchField<Type>::New
                (
                    pfld,
                    basePatch,
                    DimensionedField<Type, surfaceMesh>::null(),
                    dummyMapper
                )
            );
        }
    }

    // Add some empty patches on remaining patches (tbd.probably processor
    // patches)
    forAll(basePatchFields, patchI)
    {
        if (patchI >= patchFields.size() || !patchFields.set(patchI))
        {
            basePatchFields.set
            (
                patchI,
                fvsPatchField<Type>::New
                (
                    emptyFvsPatchField<Type>::typeName,
                    baseMesh_.boundary()[patchI],
                    DimensionedField<Type, surfaceMesh>::null()
                )
            );
        }
    }

    // Construct a volField
    IOobject baseIO
    (
        fld.name(),
        baseMesh_.time().timeName(),
        fld.local(),
        baseMesh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    );

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tfield
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            baseIO,
            baseMesh_,
            fld.dimensions(),
            internalField,
            basePatchFields
        )
    );

    tfield.ref().oriented() = fld.oriented();

    return tfield;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::parFvFieldReconstructor::reconstructFvSurfaceField
(
    const IOobject& fieldIoObject
) const
{
    // Read the field
    GeometricField<Type, fvsPatchField, surfaceMesh> fld
    (
        fieldIoObject,
        procMesh_
    );

    return reconstructFvSurfaceField(fld);
}


template<class Type>
void Foam::parFvFieldReconstructor::reconstructFvVolumeInternalFields
(
    const IOobjectList& objects,
    const HashSet<word>& selectedFields
) const
{
    const word& fieldClassName = DimensionedField<Type, volMesh>::typeName;

    IOobjectList fields = objects.lookupClass(fieldClassName);

    if (fields.size())
    {
        Info<< "    Reconstructing " << fieldClassName << "s\n" << endl;

        forAllConstIter(IOobjectList, fields, fieldIter)
        {
            if
            (
                selectedFields.empty()
             || selectedFields.found(fieldIter()->name())
            )
            {
                Info<< "        " << fieldIter()->name() << endl;

                tmp<DimensionedField<Type, volMesh>> tfld
                (
                    reconstructFvVolumeInternalField<Type>(*fieldIter())
                );

                if (isWriteProc_)
                {
                    tfld().write();
                }
            }
        }
        Info<< endl;
    }
}


template<class Type>
void Foam::parFvFieldReconstructor::reconstructFvVolumeFields
(
    const IOobjectList& objects,
    const HashSet<word>& selectedFields
) const
{
    const word& fieldClassName =
        GeometricField<Type, fvPatchField, volMesh>::typeName;

    IOobjectList fields = objects.lookupClass(fieldClassName);

    if (fields.size())
    {
        Info<< "    Reconstructing " << fieldClassName << "s\n" << endl;

        forAllConstIter(IOobjectList, fields, fieldIter)
        {
            const word& name = fieldIter()->name();

            if
            (
                (selectedFields.empty() || selectedFields.found(name))
             && name != "cellDist"
            )
            {
                Info<< "        " << name << endl;

                tmp<GeometricField<Type, fvPatchField, volMesh>> tfld
                (
                    reconstructFvVolumeField<Type>(*fieldIter())
                );
                if (isWriteProc_)
                {
                    tfld().write();
                }
            }
        }
        Info<< endl;
    }
}


template<class Type>
void Foam::parFvFieldReconstructor::reconstructFvSurfaceFields
(
    const IOobjectList& objects,
    const HashSet<word>& selectedFields
) const
{
    const word& fieldClassName =
        GeometricField<Type, fvsPatchField, surfaceMesh>::typeName;

    IOobjectList fields = objects.lookupClass(fieldClassName);

    if (fields.size())
    {
        Info<< "    Reconstructing " << fieldClassName << "s\n" << endl;

        forAllConstIter(IOobjectList, fields, fieldIter)
        {
            if
            (
                selectedFields.empty()
             || selectedFields.found(fieldIter()->name())
            )
            {
                Info<< "        " << fieldIter()->name() << endl;

                tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tfld
                (
                    reconstructFvSurfaceField<Type>(*fieldIter())
                );
                if (isWriteProc_)
                {
                    tfld().write();
                }
            }
        }
        Info<< endl;
    }
}


// ************************************************************************* //
