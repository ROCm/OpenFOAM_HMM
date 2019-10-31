/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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

#include "pointFieldReconstructor.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::pointPatchField, Foam::pointMesh>>
Foam::pointFieldReconstructor::reconstructField(const IOobject& fieldIoObject)
{
    // Read the field for all the processors
    PtrList<GeometricField<Type, pointPatchField, pointMesh>> procFields
    (
        procMeshes_.size()
    );

    forAll(procMeshes_, proci)
    {
        procFields.set
        (
            proci,
            new GeometricField<Type, pointPatchField, pointMesh>
            (
                IOobject
                (
                    fieldIoObject.name(),
                    procMeshes_[proci]().time().timeName(),
                    procMeshes_[proci](),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                procMeshes_[proci]
            )
        );
    }


    // Create the internalField
    Field<Type> internalField(mesh_.size());

    // Create the patch fields
    PtrList<pointPatchField<Type>> patchFields(mesh_.boundary().size());


    forAll(procMeshes_, proci)
    {
        const GeometricField<Type, pointPatchField, pointMesh>&
            procField = procFields[proci];

        // Get processor-to-global addressing for use in rmap
        const labelList& procToGlobalAddr = pointProcAddressing_[proci];

        // Set the cell values in the reconstructed field
        internalField.rmap
        (
            procField.primitiveField(),
            procToGlobalAddr
        );

        // Set the boundary patch values in the reconstructed field
        forAll(boundaryProcAddressing_[proci], patchi)
        {
            // Get patch index of the original patch
            const label curBPatch = boundaryProcAddressing_[proci][patchi];

            // check if the boundary patch is not a processor patch
            if (curBPatch >= 0)
            {
                if (!patchFields(curBPatch))
                {
                    patchFields.set(
                        curBPatch,
                        pointPatchField<Type>::New
                        (
                            procField.boundaryField()[patchi],
                            mesh_.boundary()[curBPatch],
                            DimensionedField<Type, pointMesh>::null(),
                            pointPatchFieldReconstructor
                            (
                                mesh_.boundary()[curBPatch].size()
                            )
                        )
                    );
                }

                patchFields[curBPatch].rmap
                (
                    procField.boundaryField()[patchi],
                    patchPointAddressing_[proci][patchi]
                );
            }
        }
    }

    // Construct and write the field
    // setting the internalField and patchFields
    return tmp<GeometricField<Type, pointPatchField, pointMesh>>::New
    (
        IOobject
        (
            fieldIoObject.name(),
            mesh_().time().timeName(),
            mesh_(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        procFields[0].dimensions(),
        internalField,
        patchFields
    );
}


template<class Type>
Foam::label Foam::pointFieldReconstructor::reconstructFields
(
    const IOobjectList& objects,
    const UList<word>& fieldNames
)
{
    typedef GeometricField<Type, pointPatchField, pointMesh> fieldType;

    label nFields = 0;
    for (const word& fieldName : fieldNames)
    {
        const IOobject* io = objects.cfindObject<fieldType>(fieldName);
        if (io)
        {
            if (!nFields++)
            {
                Info<< "    Reconstructing "
                    << fieldType::typeName << "s\n" << nl;
            }
            Info<< "        " << fieldName << endl;

            reconstructField<Type>(*io)().write();
            ++nReconstructed_;
        }
    }

    if (nFields) Info<< endl;
    return nFields;
}


template<class Type>
Foam::label Foam::pointFieldReconstructor::reconstructFields
(
    const IOobjectList& objects,
    const wordRes& selectedFields
)
{
    typedef GeometricField<Type, pointPatchField, pointMesh> fieldType;

    const wordList fieldNames =
    (
        selectedFields.empty()
      ? objects.sortedNames<fieldType>()
      : objects.sortedNames<fieldType>(selectedFields)
    );

    return reconstructFields<Type>(objects, fieldNames);
}


template<class Type>
Foam::label Foam::pointFieldReconstructor::reconstructFields
(
    const IOobjectList& objects,
    const wordHashSet& selectedFields
)
{
    typedef GeometricField<Type, pointPatchField, pointMesh> fieldType;

    const wordList fieldNames =
    (
        selectedFields.empty()
      ? objects.sortedNames<fieldType>()
      : objects.sortedNames<fieldType>(selectedFields)
    );

    return reconstructFields<Type>(objects, fieldNames);
}


// ************************************************************************* //
