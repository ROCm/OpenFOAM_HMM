/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
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

#include "parLagrangianRedistributor.H"
#include "Time.H"
#include "IOobjectList.H"
#include "mapDistributePolyMesh.H"
#include "cloud.H"
#include "CompactIOField.H"
#include "passivePositionParticleCloud.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Container>
Foam::wordList Foam::parLagrangianRedistributor::filterObjects
(
    const IOobjectList& objects,
    const wordRes& selectedFields
)
{
    wordList fieldNames =
    (
        selectedFields.empty()
      ? objects.names<Container>()
      : objects.names<Container>(selectedFields)
    );

    // Parallel synchronise
    // - Combine names from all processors

    Pstream::combineGather(fieldNames, ListOps::uniqueEqOp<word>());
    Pstream::combineScatter(fieldNames);

    // Sort for consistent order on all processors
    Foam::sort(fieldNames);

    return fieldNames;
}


template<class Type>
Foam::label Foam::parLagrangianRedistributor::redistributeFields
(
    const mapDistributeBase& map,
    const word& cloudName,
    const IOobjectList& objects,
    const wordRes& selectedFields
) const
{
    typedef IOField<Type> fieldType;

    const wordList fieldNames
    (
        filterObjects<IOField<Type>>
        (
            objects,
            selectedFields
        )
    );

    label nFields = 0;
    for (const word& objectName : fieldNames)
    {
        if (!nFields)
        {
            Info<< "    Redistributing lagrangian "
                << fieldType::typeName << "s\n" << endl;
        }
        Info<< "        " <<  objectName << endl;

        // Read if present
        IOField<Type> field
        (
            IOobject
            (
                objectName,
                srcMesh_.time().timeName(),
                cloud::prefix/cloudName,
                srcMesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            ),
            label(0)
        );

        map.distribute(field);


        const IOobject fieldIO
        (
            objectName,
            tgtMesh_.time().timeName(),
            cloud::prefix/cloudName,
            tgtMesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        );

        if (field.size())
        {
            IOField<Type>
            (
                fieldIO,
                std::move(field)
            ).write();
        }
        else
        {
            // When running with -overwrite it should also delete the old
            // files. Below works but is not optimal.

            const fileName fldName(fieldIO.objectPath());
            Foam::rm(fldName);
        }
    }

    return nFields;
}


template<class Type>
Foam::label Foam::parLagrangianRedistributor::redistributeFieldFields
(
    const mapDistributeBase& map,
    const word& cloudName,
    const IOobjectList& objects,
    const wordRes& selectedFields
) const
{
    typedef CompactIOField<Field<Type>, Type> fieldType;

    wordList fieldNames
    (
        filterObjects<fieldType>
        (
            objects,
            selectedFields
        )
    );

    // Append IOField Field names
    {
        wordList ioFieldNames
        (
            filterObjects<IOField<Field<Type>>>
            (
                objects,
                selectedFields
            )
        );
        fieldNames.append(ioFieldNames);
    }

    label nFields = 0;
    for (const word& objectName : fieldNames)
    {
        if (!nFields++)
        {
            Info<< "    Redistributing lagrangian "
                << fieldType::typeName << "s\n" << nl;
        }
        Info<< "        " <<  objectName << nl;

        // Read if present
        CompactIOField<Field<Type>, Type> field
        (
            IOobject
            (
                objectName,
                srcMesh_.time().timeName(),
                cloud::prefix/cloudName,
                srcMesh_,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            ),
            label(0)
        );

        // Distribute
        map.distribute(field);

        // Write
        const IOobject fieldIO
        (
            objectName,
            tgtMesh_.time().timeName(),
            cloud::prefix/cloudName,
            tgtMesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        );

        if (field.size())
        {
            CompactIOField<Field<Type>, Type>
            (
                fieldIO,
                std::move(field)
            ).write();
        }
        else
        {
            // When running with -overwrite it should also delete the old
            // files. Below works but is not optimal.

            const fileName fldName(fieldIO.objectPath());
            Foam::rm(fldName);
        }
    }

    if (nFields) Info<< endl;
    return nFields;
}


template<class Container>
Foam::label Foam::parLagrangianRedistributor::readFields
(
    const passivePositionParticleCloud& cloud,
    const IOobjectList& objects,
    const wordRes& selectedFields
)
{
    const word fieldClassName(Container::typeName);

    const wordList fieldNames
    (
        filterObjects<Container>
        (
            objects,
            selectedFields
        )
    );

    label nFields = 0;
    for (const word& objectName : fieldNames)
    {
        if (!nFields++)
        {
            Info<< "    Reading lagrangian "
                << Container::typeName << "s\n" << nl;
        }
        Info<< "        " <<  objectName << nl;

        // Read if present
        Container* fieldPtr = new Container
        (
            IOobject
            (
                objectName,
                cloud.time().timeName(),
                cloud,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            label(0)
        );

        fieldPtr->store();
    }

    return nFields;
}


template<class Container>
Foam::label Foam::parLagrangianRedistributor::redistributeStoredFields
(
    const mapDistributeBase& map,
    passivePositionParticleCloud& cloud
) const
{
    HashTable<Container*> fields
    (
        cloud.lookupClass<Container>()
    );

    label nFields = 0;
    forAllIters(fields, iter)
    {
        Container& field = *(iter.val());

        if (!nFields++)
        {
            Info<< "    Redistributing lagrangian "
                << Container::typeName << "s\n" << endl;
        }
        Info<< "        " <<  field.name() << endl;

        map.distribute(field);

        const IOobject fieldIO
        (
            field.name(),
            tgtMesh_.time().timeName(),
            cloud::prefix/cloud.name(),
            tgtMesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        );

        if (field.size())
        {
            Container
            (
                fieldIO,
                std::move(field)
            ).write();
        }
        else
        {
            // When running with -overwrite it should also delete the old
            // files. Below works but is not optimal.

            const fileName fldName(fieldIO.objectPath());
            Foam::rm(fldName);
        }
    }

    return nFields;
}


// ************************************************************************* //
