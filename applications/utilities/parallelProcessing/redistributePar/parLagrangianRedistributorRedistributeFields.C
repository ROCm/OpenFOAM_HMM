/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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
    const wordHashSet& selectedFields
)
{
    // Parallel synchronise
    wordList fieldNames =
    (
        selectedFields.empty()
      ? objects.names(Container::typeName)
      : objects.names(Container::typeName, selectedFields)
    );

    Pstream::combineGather(fieldNames, ListOps::uniqueEqOp<word>());
    Pstream::combineScatter(fieldNames);

    // Ensure order is consistent
    Foam::sort(fieldNames);

    return fieldNames;
}


template<class Type>
void Foam::parLagrangianRedistributor::redistributeLagrangianFields
(
    const mapDistributeBase& map,
    const word& cloudName,
    const IOobjectList& objects,
    const wordHashSet& selectedFields
) const
{
    const word fieldClassName(IOField<Type>::typeName);

    const wordList objectNames
    (
        filterObjects<IOField<Type>>
        (
            objects,
            selectedFields
        )
    );

    if (objectNames.size())
    {
        Info<< "    Redistributing lagrangian "
            << fieldClassName << "s\n" << endl;
    }

    for (const word& objectName : objectNames)
    {
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


        if (field.size())
        {
            IOField<Type>
            (
                IOobject
                (
                    objectName,
                    tgtMesh_.time().timeName(),
                    cloud::prefix/cloudName,
                    tgtMesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                std::move(field)
            ).write();
        }
    }

    if (objectNames.size()) Info<< endl;
}


template<class Type>
void Foam::parLagrangianRedistributor::redistributeLagrangianFieldFields
(
    const mapDistributeBase& map,
    const word& cloudName,
    const IOobjectList& objects,
    const wordHashSet& selectedFields
) const
{
    const word fieldClassName(CompactIOField<Field<Type>, Type>::typeName);

    wordList objectNames
    (
        filterObjects<CompactIOField<Field<Type>, Type>>
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
        objectNames.append(ioFieldNames);
    }


    if (objectNames.size())
    {
        Info<< "    Redistributing lagrangian "
            << fieldClassName << "s\n" << endl;
    }

    for (const word& objectName : objectNames)
    {
        Info<< "        " <<  objectName << endl;

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
        if (field.size())
        {
            CompactIOField<Field<Type>, Type>
            (
                IOobject
                (
                    objectName,
                    tgtMesh_.time().timeName(),
                    cloud::prefix/cloudName,
                    tgtMesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                std::move(field)
            ).write();
        }
    }
}


template<class Container>
void Foam::parLagrangianRedistributor::readLagrangianFields
(
    const passivePositionParticleCloud& cloud,
    const IOobjectList& objects,
    const wordHashSet& selectedFields
)
{
    const word fieldClassName(Container::typeName);

    const wordList objectNames
    (
        filterObjects<Container>
        (
            objects,
            selectedFields
        )
    );

    if (objectNames.size())
    {
        Info<< "    Reading lagrangian "
            << fieldClassName << "s\n" << endl;
    }

    for (const word& objectName : objectNames)
    {
        Info<< "        " <<  objectName << endl;

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
}


template<class Container>
void Foam::parLagrangianRedistributor::redistributeStoredLagrangianFields
(
    const mapDistributeBase& map,
    passivePositionParticleCloud& cloud
) const
{
    const word fieldClassName(Container::typeName);

    HashTable<Container*> fields
    (
        cloud.lookupClass<Container>()
    );

    if (fields.size())
    {
        Info<< "    Redistributing lagrangian "
            << fieldClassName << "s\n" << endl;
    }

    forAllIters(fields, iter)
    {
        Container& field = *(*iter);

        Info<< "        " <<  field.name() << endl;

        map.distribute(field);

        if (field.size())
        {
            Container
            (
                IOobject
                (
                    field.name(),
                    tgtMesh_.time().timeName(),
                    cloud::prefix/cloud.name(),
                    tgtMesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                std::move(field)
            ).write();
        }
    }
}


// ************************************************************************* //
