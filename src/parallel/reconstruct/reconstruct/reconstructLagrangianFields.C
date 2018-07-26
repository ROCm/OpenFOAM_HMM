/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2018 OpenCFD Ltd.
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

#include "IOField.H"
#include "CompactIOField.H"
#include "Time.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::IOField<Type>> Foam::reconstructLagrangianField
(
    const word& cloudName,
    const polyMesh& mesh,
    const PtrList<fvMesh>& meshes,
    const word& fieldName
)
{
    // Construct empty field on mesh
    auto tfield = tmp<IOField<Type>>::New
    (
        IOobject
        (
            fieldName,
            mesh.time().timeName(),
            cloud::prefix/cloudName,
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Field<Type>(0)
    );
    auto& field = tfield.ref();

    for (const fvMesh& localMesh : meshes)
    {
        // Check object on local mesh
        IOobject localIOobject
        (
            fieldName,
            localMesh.time().timeName(),
            cloud::prefix/cloudName,
            localMesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if (localIOobject.typeHeaderOk<IOField<Type>>(true))
        {
            IOField<Type> localField(localIOobject);

            const label offset = field.size();
            field.setSize(offset + localField.size());

            forAll(localField, j)
            {
                field[offset + j] = localField[j];
            }
        }
    }

    return tfield;
}


template<class Type>
Foam::tmp<Foam::CompactIOField<Foam::Field<Type>, Type>>
Foam::reconstructLagrangianFieldField
(
    const word& cloudName,
    const polyMesh& mesh,
    const PtrList<fvMesh>& meshes,
    const word& fieldName
)
{
    // Construct empty field on mesh
    auto tfield = tmp<CompactIOField<Field<Type>, Type>>::New
    (
        IOobject
        (
            fieldName,
            mesh.time().timeName(),
            cloud::prefix/cloudName,
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Field<Field<Type>>(0)
    );
    auto& field = tfield.ref();

    for (const fvMesh& localMesh : meshes)
    {
        // Check object on local mesh
        IOobject localIOobject
        (
            fieldName,
            localMesh.time().timeName(),
            cloud::prefix/cloudName,
            localMesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if
        (
            localIOobject.typeHeaderOk<CompactIOField<Field<Type>, Type>>
            (
                false
            )
         || localIOobject.typeHeaderOk<IOField<Field<Type>>>(false)
        )
        {
            CompactIOField<Field<Type>, Type> localField(localIOobject);

            const label offset = field.size();
            field.setSize(offset + localField.size());

            forAll(localField, j)
            {
                field[offset + j] = localField[j];
            }
        }
    }

    return tfield;
}


template<class Type>
void Foam::reconstructLagrangianFields
(
    const word& cloudName,
    const polyMesh& mesh,
    const PtrList<fvMesh>& meshes,
    const IOobjectList& objects,
    const wordHashSet& selectedFields
)
{
    const word& clsName = IOField<Type>::typeName;

    const wordList fieldNames =
    (
        selectedFields.empty()
      ? objects.sortedNames(clsName)
      : objects.sortedNames(clsName, selectedFields)
    );

    if (fieldNames.size())
    {
        Info<< "    Reconstructing lagrangian " << clsName << "s\n" << nl;
    }

    for (const word& fieldName : fieldNames)
    {
        Info<< "        " << fieldName << endl;

        reconstructLagrangianField<Type>
        (
            cloudName,
            mesh,
            meshes,
            fieldName
        )().write();
    }

    if (fieldNames.size()) Info<< endl;
}


template<class Type>
void Foam::reconstructLagrangianFieldFields
(
    const word& cloudName,
    const polyMesh& mesh,
    const PtrList<fvMesh>& meshes,
    const IOobjectList& objects,
    const wordHashSet& selectedFields
)
{
    {
        const word& clsName = CompactIOField<Field<Type>,Type>::typeName;

        const wordList fieldNames =
        (
            selectedFields.empty()
          ? objects.sortedNames(clsName)
          : objects.sortedNames(clsName, selectedFields)
        );

        if (fieldNames.size())
        {
            Info<< "    Reconstructing lagrangian " << clsName << "s\n" << nl;
        }

        for (const word& fieldName : fieldNames)
        {
            Info<< "        " << fieldName << endl;

            reconstructLagrangianFieldField<Type>
            (
                cloudName,
                mesh,
                meshes,
                fieldName
            )().write();
        }

        if (fieldNames.size()) Info<< endl;
    }

    {
        const word& clsName = IOField<Field<Type>>::typeName;

        const wordList fieldNames =
        (
            selectedFields.empty()
          ? objects.sortedNames(clsName)
          : objects.sortedNames(clsName, selectedFields)
        );

        if (fieldNames.size())
        {
            Info<< "    Reconstructing lagrangian " << clsName << "s\n" << nl;
        }

        for (const word& fieldName : fieldNames)
        {
            Info<< "        " << fieldName << endl;

            reconstructLagrangianFieldField<Type>
            (
                cloudName,
                mesh,
                meshes,
                fieldName
            )().write();
        }

        if (fieldNames.size()) Info<< endl;
    }
}


// ************************************************************************* //
