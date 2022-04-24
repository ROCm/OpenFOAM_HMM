/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class GeoField>
void Foam::fieldsDistributor::readField
(
    const IOobject& io,
    const typename GeoField::Mesh& mesh,
    const label i,
    PtrList<GeoField>& fields
)
{
    fields.set(i, new GeoField(io, mesh));
}

template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::fieldsDistributor::readField
(
    const IOobject& io,
    const typename GeoMesh::Mesh& mesh,
    const label i,
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& fields
)
{
    fields.set
    (
        i,
        new GeometricField<Type, PatchField, GeoMesh>(io, mesh, false)
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
void Foam::fieldsDistributor::readFields
(
    const typename GeoMesh::Mesh& mesh,
    const IOobjectList& objects,
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& fields,
    const bool readOldTime
)
{
    typedef GeometricField<Type, PatchField, GeoMesh> GeoField;

    // GeoField fields - sorted for consistent order on all processors
    UPtrList<const IOobject> fieldObjects(objects.sorted<GeoField>());

    // Construct the fields
    fields.resize(fieldObjects.size());

    forAll(fieldObjects, i)
    {
        fields.set(i, new GeoField(fieldObjects[i], mesh, readOldTime));
    }
}


template<class Mesh, class GeoField>
void Foam::fieldsDistributor::readFields
(
    const Mesh& mesh,
    const IOobjectList& objects,
    PtrList<GeoField>& fields
)
{
    // GeoField fields - sorted for consistent order on all processors
    UPtrList<const IOobject> fieldObjects(objects.sorted<GeoField>());

    // Construct the fields
    fields.resize(fieldObjects.size());

    forAll(fieldObjects, i)
    {
        fields.set(i, new GeoField(fieldObjects[i], mesh));
    }
}


template<class GeoField, class MeshSubsetter>
void Foam::fieldsDistributor::readFields
(
    const boolList& haveMeshOnProc,
    const typename GeoField::Mesh& mesh,
    const autoPtr<MeshSubsetter>& subsetterPtr,
    IOobjectList& allObjects,
    PtrList<GeoField>& fields,
    const bool deregister
)
{
    // Get my objects of type
    IOobjectList objects(allObjects.lookupClass<GeoField>());

    // Check that we all have all objects
    wordList objectNames = objects.sortedNames();

    // Get master names
    wordList masterNames(objectNames);
    Pstream::broadcast(masterNames);

    if (haveMeshOnProc[Pstream::myProcNo()] && objectNames != masterNames)
    {
        FatalErrorInFunction
            << "Objects not synchronised across processors." << nl
            << "Master has " << flatOutput(masterNames) << nl
            << "Processor " << Pstream::myProcNo()
            << " has " << flatOutput(objectNames)
            << exit(FatalError);
    }

    fields.clear();
    fields.resize(masterNames.size());

    if (fields.empty())
    {
        if (deregister)
        {
            // Extra safety - remove all such types
            HashTable<const GeoField*> other
            (
                mesh.thisDb().objectRegistry::template lookupClass<GeoField>()
            );

            forAllConstIters(other, iter)
            {
                GeoField& fld = const_cast<GeoField&>(*iter.val());

                if (!fld.ownedByRegistry())
                {
                    fld.checkOut();
                }
            }
        }

        // Early exit
        return;
    }


    // Have master send all fields to processors that don't have a mesh. The
    // issue is if a patchField does any parallel operations inside its
    // construct-from-dictionary. This will not work when going to more
    // processors (e.g. decompose = 1 -> many) ! We could make a special
    // exception for decomposePar but nicer would be to have read-communicator
    // ... For now detect if decomposing & disable parRun
    if (Pstream::master())
    {
        // Work out if we're decomposing - none of the subprocs has a mesh
        bool decompose = true;
        for (const int proci : Pstream::subProcs())
        {
            if (haveMeshOnProc[proci])
            {
                decompose = false;
            }
        }

        const bool oldParRun = Pstream::parRun();
        if (decompose)
        {
            Pstream::parRun(false);
        }

        forAll(masterNames, i)
        {
            const word& name = masterNames[i];
            IOobject& io = *objects[name];
            io.writeOpt(IOobject::AUTO_WRITE);

            // Load field (but not oldTime)
            readField(io, mesh, i, fields);
        }

        Pstream::parRun(oldParRun);  // Restore any changes
    }
    else if (haveMeshOnProc[Pstream::myProcNo()])
    {
        // Have mesh so just try to load
        forAll(masterNames, i)
        {
            const word& name = masterNames[i];
            IOobject& io = *objects[name];
            io.writeOpt(IOobject::AUTO_WRITE);

            /// Pout<< "Attempt read: " << name << endl;

            // Load field (but not oldTime)
            readField(io, mesh, i, fields);
        }
    }


    // Missing fields on any processors?
    // - construct from dictionary

    PtrList<dictionary> fieldDicts;

    if (Pstream::master())
    {
        // Broadcast zero sized fields everywhere (if needed)
        // Send like a list of dictionaries

        OPBstream toProcs(UPstream::masterNo());  // worldComm

        const label nDicts = (subsetterPtr ? fields.size() : label(0));

        toProcs << nDicts << token::BEGIN_LIST;  // Begin list

        if (nDicts)
        {
            // Disable communication for interpolate() method
            const bool oldParRun = Pstream::parRun(false);

            const auto& subsetter = subsetterPtr();

            forAll(fields, i)
            {
                tmp<GeoField> tsubfld = subsetter.interpolate(fields[i]);

                // Surround each with {} as dictionary entry
                toProcs.beginBlock();
                toProcs << tsubfld();
                toProcs.endBlock();
            }

            Pstream::parRun(oldParRun);  // Restore state
        }

        toProcs << token::END_LIST << token::NL;  // End list
    }
    else
    {
        // Receive the broadcast...
        IPBstream fromMaster(UPstream::masterNo());  // worldComm

        // But only consume where needed...
        if (!haveMeshOnProc[Pstream::myProcNo()])
        {
            fromMaster >> fieldDicts;
        }
    }


    // Use the received dictionaries to create fields
    // (will be empty if we didn't require them)

    // Disable communication when constructing from dictionary
    const bool oldParRun = Pstream::parRun(false);

    forAll(fieldDicts, i)
    {
        fields.set
        (
            i,
            new GeoField
            (
                IOobject
                (
                    masterNames[i],
                    mesh.time().timeName(),
                    mesh.thisDb(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                fieldDicts[i]
            )
        );
    }

    Pstream::parRun(oldParRun);  // Restore any changes


    // Finally. Can checkOut of registry as required
    if (deregister)
    {
        /// Info<< "De-registering fields:";
        for (auto& fld : fields)
        {
            /// Info<< "  " << fld.name();

            // Ensure it is not destroyed by polyMesh deletion
            fld.checkOut();
        }
        /// Info<< nl;

        // Extra safety - remove all such types
        HashTable<const GeoField*> other
        (
            mesh.thisDb().objectRegistry::template lookupClass<GeoField>()
        );

        forAllConstIters(other, iter)
        {
            GeoField& fld = const_cast<GeoField&>(*iter.val());

            if (!fld.ownedByRegistry())
            {
                fld.checkOut();
            }
        }
    }
}


// ************************************************************************* //
