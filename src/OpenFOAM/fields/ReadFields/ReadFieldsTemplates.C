/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "ReadFields.H"
#include "HashSet.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
Foam::wordList Foam::ReadFields
(
    const typename GeoMesh::Mesh& mesh,
    const IOobjectList& objects,
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& fields,
    const bool syncPar,
    const bool readOldTime
)
{
    typedef GeometricField<Type, PatchField, GeoMesh> GeoField;

    // Names of GeoField objects, sorted order.
    const wordList fieldNames(objects.names(GeoField::typeName, syncPar));

    // Construct the fields - reading in consistent (master) order.
    fields.resize(fieldNames.size());

    label nFields = 0;

    for (const word& fieldName : fieldNames)
    {
        if (!nFields)
        {
            Info<< "Reading " << GeoField::typeName << ':';
        }
        Info<< ' ' << fieldName;

        const IOobject& io = *objects[fieldName];

        fields.set
        (
            nFields++,
            new GeoField
            (
                IOobject
                (
                    io.name(),
                    io.instance(),
                    io.local(),
                    io.db(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE,
                    io.registerObject()
                ),
                mesh,
                readOldTime
            )
        );
    }

    if (nFields) Info<< endl;

    return fieldNames;
}


template<class GeoField, class Mesh>
Foam::wordList Foam::ReadFields
(
    const Mesh& mesh,
    const IOobjectList& objects,
    PtrList<GeoField>& fields,
    const bool syncPar
)
{
    // Names of GeoField objects, sorted order.
    const wordList fieldNames(objects.names(GeoField::typeName, syncPar));

    // Construct the fields - reading in consistent (master) order.
    fields.resize(fieldNames.size());

    label nFields = 0;

    for (const word& fieldName : fieldNames)
    {
        if (!nFields)
        {
            Info<< "Reading " << GeoField::typeName << ':';
        }
        Info<< ' ' << fieldName;

        const IOobject& io = *objects[fieldName];

        fields.set
        (
            nFields++,
            new GeoField
            (
                IOobject
                (
                    io.name(),
                    io.instance(),
                    io.local(),
                    io.db(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE,
                    io.registerObject()
                ),
                mesh
            )
        );
    }

    if (nFields) Info<< endl;

    return fieldNames;
}


template<class GeoField>
Foam::wordList Foam::ReadFields
(
    const IOobjectList& objects,
    PtrList<GeoField>& fields,
    const bool syncPar
)
{
    // Names of GeoField objects, sorted order.
    const wordList fieldNames(objects.names(GeoField::typeName, syncPar));

    // Construct the fields - reading in consistent (master) order.
    fields.resize(fieldNames.size());

    label nFields = 0;

    for (const word& fieldName : fieldNames)
    {
        if (!nFields)
        {
            Info<< "Reading " << GeoField::typeName << ':';
        }
        Info<< ' ' << fieldName;

        const IOobject& io = *objects[fieldName];

        fields.set
        (
            nFields++,
            new GeoField
            (
                IOobject
                (
                    io.name(),
                    io.instance(),
                    io.local(),
                    io.db(),
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE,
                    io.registerObject()
                )
            )
        );
    }

    if (nFields) Info<< endl;

    return fieldNames;
}


template<class GeoField>
void Foam::ReadFields
(
    const word& fieldName,
    const typename GeoField::Mesh& mesh,
    const wordList& timeNames,
    objectRegistry& fieldsCache
)
{
    // Unload times that are no longer used
    {
        wordHashSet unusedTimes(fieldsCache.toc());
        unusedTimes.erase(timeNames);

        //Info<< "Unloading times " << unusedTimes << endl;

        for (const word& timeName : unusedTimes)
        {
            objectRegistry& timeCache =
                fieldsCache.lookupObjectRef<objectRegistry>(timeName);

            fieldsCache.checkOut(timeCache);
        }
    }


    // Load any new fields
    for (const word& timeName : timeNames)
    {
        // Create if not found
        if (!fieldsCache.found(timeName))
        {
            //Info<< "Creating registry for time " << timeName << endl;

            // Create objectRegistry if not found
            objectRegistry* timeCachePtr = new objectRegistry
            (
                IOobject
                (
                    timeName,
                    timeName,
                    fieldsCache,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                )
            );
            timeCachePtr->store();
        }

        // Obtain cache for current time
        const objectRegistry& timeCache =
            fieldsCache.lookupObject<objectRegistry>(timeName);

        // Store field if not found
        if (!timeCache.found(fieldName))
        {
            //Info<< "Loading field " << fieldName
            //    << " for time " << timeName << endl;

            GeoField loadedFld
            (
                IOobject
                (
                    fieldName,
                    timeName,
                    mesh.thisDb(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false  // do not register
                ),
                mesh
            );

            // Transfer to timeCache (new objectRegistry and store flag)
            GeoField* fldPtr = new GeoField
            (
                IOobject
                (
                    fieldName,
                    timeName,
                    timeCache,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                loadedFld
            );
            fldPtr->store();
        }
    }
}


template<class GeoField>
void Foam::ReadFields
(
    const word& fieldName,
    const typename GeoField::Mesh& mesh,
    const wordList& timeNames,
    const word& registryName
)
{
    ReadFields<GeoField>
    (
        fieldName,
        mesh,
        timeNames,
        const_cast<objectRegistry&>
        (
            mesh.thisDb().subRegistry(registryName, true)
        )
    );
}


template<class GeoFieldType>
void Foam::readFields
(
    const typename GeoFieldType::Mesh& mesh,
    const IOobjectList& objects,
    const wordHashSet& selectedFields,
    LIFOStack<regIOobject*>& storedObjects
)
{
    // Names of GeoField objects, sorted order. Not synchronised.
    const wordList fieldNames
    (
        objects.sortedNames
        (
            GeoFieldType::typeName,
            selectedFields   // Only permit these
        )
    );

    label nFields = 0;

    for (const word& fieldName : fieldNames)
    {
        const IOobject& io = *objects[fieldName];

        if (!nFields)
        {
            Info<< "    " << GeoFieldType::typeName << ':';
        }
        Info<< ' ' << fieldName;

        GeoFieldType* fieldPtr = new GeoFieldType
        (
            IOobject
            (
                fieldName,
                io.instance(),
                io.local(),
                io.db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
        fieldPtr->store();
        storedObjects.push(fieldPtr);

        ++nFields;
    }

    if (nFields) Info<< endl;
}


template<class UniformFieldType>
void Foam::readUniformFields
(
    const IOobjectList& objects,
    const wordHashSet& selectedFields,
    LIFOStack<regIOobject*>& storedObjects,
    const bool syncPar
)
{
    // Names of UniformField objects, sorted order.
    const wordList fieldNames
    (
        objects.names
        (
            UniformFieldType::typeName,
            selectedFields,  // Only permit these
            syncPar
        )
    );

    label nFields = 0;

    for (const word& fieldName : fieldNames)
    {
        const IOobject& io = *objects[fieldName];

        if (!nFields)
        {
            Info<< "    " << UniformFieldType::typeName << ':';
        }
        Info<< ' ' << fieldName;

        UniformFieldType* fieldPtr = new UniformFieldType
        (
            IOobject
            (
                fieldName,
                io.instance(),
                io.local(),
                io.db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );
        fieldPtr->store();
        storedObjects.push(fieldPtr);

        ++nFields;
    }

    if (nFields) Info<< endl;
}


// ************************************************************************* //
