/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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

#include "sampledSurfaces.H"
#include "polySurface.H"

#include "mapPolyMesh.H"
#include "volFields.H"
#include "HashOps.H"
#include "ListOps.H"
#include "Time.H"
#include "UIndirectList.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledSurfaces, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        sampledSurfaces,
        dictionary
    );
}

Foam::scalar Foam::sampledSurfaces::mergeTol_ = 1e-10;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::polySurface* Foam::sampledSurfaces::getRegistrySurface
(
    const sampledSurface& s
) const
{
    return s.getRegistrySurface
    (
        storedObjects(),
        IOobject::groupName(name(), s.name())
    );
}


Foam::polySurface* Foam::sampledSurfaces::storeRegistrySurface
(
    const sampledSurface& s
)
{
    return s.storeRegistrySurface
    (
        storedObjects(),
        IOobject::groupName(name(), s.name())
    );
}


bool Foam::sampledSurfaces::removeRegistrySurface
(
    const sampledSurface& s
)
{
    return s.removeRegistrySurface
    (
        storedObjects(),
        IOobject::groupName(name(), s.name())
    );
}


void Foam::sampledSurfaces::countFields()
{
    wordList allFields;    // Just needed for warnings
    HashTable<wordHashSet> selected;

    if (loadFromFiles_)
    {
        // Check files for a particular time
        IOobjectList objects(obr_, obr_.time().timeName());

        allFields = objects.names();
        selected = objects.classes(fieldSelection_);
    }
    else
    {
        // Check currently available fields
        allFields = obr_.names();
        selected = obr_.classes(fieldSelection_);
    }

    if (Pstream::parRun())
    {
        Pstream::mapCombineGather(selected, HashSetOps::plusEqOp<word>());
        Pstream::mapCombineScatter(selected);
    }


    DynamicList<label> missed(fieldSelection_.size());

    // Detect missing fields
    forAll(fieldSelection_, i)
    {
        if (findStrings(fieldSelection_[i], allFields).empty())
        {
            missed.append(i);
        }
    }

    if (missed.size())
    {
        WarningInFunction
            << nl
            << "Cannot find "
            << (loadFromFiles_ ? "field file" : "registered field")
            << " matching "
            << UIndirectList<wordRe>(fieldSelection_, missed) << endl;
    }


    // Currently only support volume and surface field types
    label nVolumeFields = 0;
    label nSurfaceFields = 0;

    forAllConstIters(selected, iter)
    {
        const word& clsName = iter.key();
        const label n = iter.val().size();

        if (fieldTypes::volume.found(clsName))
        {
            nVolumeFields += n;
        }
        else if (sampledSurface::surfaceFieldTypes.found(clsName))
        {
            nSurfaceFields += n;
        }
    }

    // Now propagate field counts (per surface)
    // - can update writer even when not writing without problem

    forAll(*this, surfi)
    {
        const sampledSurface& s = (*this)[surfi];
        surfaceWriter& outWriter = writers_[surfi];

        outWriter.nFields
        (
            nVolumeFields
          + (s.withSurfaceFields() ? nSurfaceFields : 0)
          +
            (
                // Face-id field, but not if the writer manages that itself
                !s.isPointData() && s.hasFaceIds() && !outWriter.usesFaceIds()
              ? 1 : 0
            )
        );
    }
}


Foam::autoPtr<Foam::surfaceWriter> Foam::sampledSurfaces::newWriter
(
    word writeType,
    const dictionary& formatOptions,
    const dictionary& surfDict
)
{
    // Per-surface adjustment
    surfDict.readIfPresent<word>("surfaceFormat", writeType);

    dictionary options = formatOptions.subOrEmptyDict(writeType);

    options.merge
    (
        surfDict.subOrEmptyDict("formatOptions").subOrEmptyDict(writeType)
    );

    return surfaceWriter::New(writeType, options);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledSurfaces::sampledSurfaces
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObjects::fvMeshFunctionObject(name, runTime, dict),
    PtrList<sampledSurface>(),
    loadFromFiles_(false),
    verbose_(false),
    onExecute_(false),
    outputPath_
    (
        time_.globalPath()/functionObject::outputPrefix/name
    ),
    fieldSelection_(),
    sampleFaceScheme_(),
    sampleNodeScheme_(),
    writers_(),
    actions_(),
    nFaces_()
{
    outputPath_.clean();  // Remove unneeded ".."

    read(dict);
}


Foam::sampledSurfaces::sampledSurfaces
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    functionObjects::fvMeshFunctionObject(name, obr, dict),
    PtrList<sampledSurface>(),
    loadFromFiles_(loadFromFiles),
    verbose_(false),
    onExecute_(false),
    outputPath_
    (
        time_.globalPath()/functionObject::outputPrefix/name
    ),
    fieldSelection_(),
    sampleFaceScheme_(),
    sampleNodeScheme_(),
    writers_(),
    actions_(),
    nFaces_()
{
    outputPath_.clean();  // Remove unneeded ".."

    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sampledSurfaces::verbose(const bool on)
{
    bool old(verbose_);
    verbose_ = on;
    return old;
}


bool Foam::sampledSurfaces::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    PtrList<sampledSurface>::clear();
    writers_.clear();
    actions_.clear();
    nFaces_.clear();
    fieldSelection_.clear();

    verbose_ = dict.getOrDefault("verbose", false);
    onExecute_ = dict.getOrDefault("sampleOnExecute", false);

    sampleFaceScheme_ =
        dict.getOrDefault<word>("sampleScheme", "cell");

    sampleNodeScheme_ =
        dict.getOrDefault<word>("interpolationScheme", "cellPoint");

    const entry* eptr = dict.findEntry("surfaces");

    // Surface writer type and format options
    const word writerType =
        (eptr ? dict.get<word>("surfaceFormat") : word::null);

    const dictionary formatOptions(dict.subOrEmptyDict("formatOptions"));

    // Store on registry?
    const bool dfltStore = dict.getOrDefault("store", false);

    if (eptr && eptr->isDict())
    {
        PtrList<sampledSurface> surfs(eptr->dict().size());

        actions_.resize(surfs.size(), ACTION_WRITE); // Default action
        writers_.resize(surfs.size());
        nFaces_.resize(surfs.size(), Zero);

        label surfi = 0;

        for (const entry& dEntry : eptr->dict())
        {
            if (!dEntry.isDict())
            {
                continue;
            }

            const dictionary& surfDict = dEntry.dict();

            autoPtr<sampledSurface> surf =
                sampledSurface::New
                (
                    dEntry.keyword(),
                    mesh_,
                    surfDict
                );

            if (!surf || !surf->enabled())
            {
                continue;
            }

            // Define the surface
            surfs.set(surfi, surf);

            // Define additional action(s)
            if (surfDict.getOrDefault("store", dfltStore))
            {
                actions_[surfi] |= ACTION_STORE;
            }
            if (surfDict.getOrDefault("surfMeshStore", false))
            {
                actions_[surfi] |= ACTION_SURF_MESH;
            }

            // Define surface writer, but do NOT yet attach a surface
            writers_.set
            (
                surfi,
                newWriter(writerType, formatOptions, surfDict)
            );

            writers_[surfi].isPointData(surfs[surfi].isPointData());

            // Use outputDir/TIME/surface-name
            writers_[surfi].useTimeDir(true);
            writers_[surfi].verbose(verbose_);

            ++surfi;
        }

        surfs.resize(surfi);
        actions_.resize(surfi);
        writers_.resize(surfi);
        surfaces().transfer(surfs);
    }
    else if (eptr)
    {
        // This is slightly trickier.
        // We want access to the individual dictionaries used for construction

        DynamicList<dictionary> capture;

        PtrList<sampledSurface> input
        (
            eptr->stream(),
            sampledSurface::iNewCapture(mesh_, capture)
        );

        PtrList<sampledSurface> surfs(input.size());

        actions_.resize(surfs.size(), ACTION_WRITE); // Default action
        writers_.resize(surfs.size());
        nFaces_.resize(surfs.size(), Zero);

        label surfi = 0;

        forAll(input, inputi)
        {
            const dictionary& surfDict = capture[inputi];

            autoPtr<sampledSurface> surf = input.release(inputi);

            if (!surf || !surf->enabled())
            {
                continue;
            }

            // Define the surface
            surfs.set(surfi, surf);

            // Define additional action(s)
            if (surfDict.getOrDefault("store", dfltStore))
            {
                actions_[surfi] |= ACTION_STORE;
            }
            if (surfDict.getOrDefault("surfMeshStore", false))
            {
                actions_[surfi] |= ACTION_SURF_MESH;
            }

            // Define surface writer, but do NOT yet attach a surface
            writers_.set
            (
                surfi,
                newWriter(writerType, formatOptions, surfDict)
            );

            writers_[surfi].isPointData(surfs[surfi].isPointData());

            // Use outputDir/TIME/surface-name
            writers_[surfi].useTimeDir(true);
            writers_[surfi].verbose(verbose_);

            ++surfi;
        }

        surfs.resize(surfi);
        actions_.resize(surfi);
        writers_.resize(surfi);
        surfaces().transfer(surfs);
    }


    const auto& surfs = surfaces();

    // Have some surfaces, so sort out which fields are needed and report

    if (surfs.size())
    {
        nFaces_.resize(surfs.size(), Zero);

        dict.readEntry("fields", fieldSelection_);
        fieldSelection_.uniq();

        forAll(*this, surfi)
        {
            const sampledSurface& s = (*this)[surfi];

            if (!surfi)
            {
                Info<< "Sampled surface:" << nl;
            }

            Info<< "    " << s.name() << " -> " << writers_[surfi].type();
            if (actions_[surfi] & ACTION_STORE)
            {
                Info<< ", store on registry ("
                    << IOobject::groupName(name(), s.name()) << ')';
            }
            if (actions_[surfi] & ACTION_SURF_MESH)
            {
                Info<< ", store as surfMesh (DEPRECATED)";
            }
            Info<< nl;
            Info<< "        ";
            s.print(Info, 0);
            Info<< nl;
        }
        Info<< nl;
    }

    if (debug && Pstream::master())
    {
        Pout<< "sample fields:" << fieldSelection_ << nl
            << "sample surfaces:" << nl << '(' << nl;

        for (const sampledSurface& s : surfaces())
        {
            Pout<< "  " << s << nl;
        }
        Pout<< ')' << endl;
    }

    // Ensure all surfaces and merge information are expired
    expire(true);

    return true;
}


bool Foam::sampledSurfaces::performAction(unsigned request)
{
    // Update surfaces and store
    bool ok = false;

    forAll(*this, surfi)
    {
        sampledSurface& s = (*this)[surfi];

        if (request & actions_[surfi])
        {
            if (s.update())
            {
                writers_[surfi].expire();
            }

            nFaces_[surfi] = returnReduce(s.faces().size(), sumOp<label>());

            ok = ok || nFaces_[surfi];


            // Store surfaces (even empty ones) otherwise we miss geometry
            // updates.
            // Any associated fields will be removed if the size changes

            if ((request & actions_[surfi]) & ACTION_STORE)
            {
                storeRegistrySurface(s);
            }

            if ((request & actions_[surfi]) & ACTION_SURF_MESH)
            {
                s.storeSurfMesh();
            }
        }
    }

    if (!ok)
    {
        // No surface with an applicable action or with faces to sample
        return true;
    }


    // Determine the per-surface number of fields, including Ids etc.
    // Only seems to be needed for VTK legacy
    countFields();


    // Update writers

    forAll(*this, surfi)
    {
        const sampledSurface& s = (*this)[surfi];

        if (((request & actions_[surfi]) & ACTION_WRITE) && nFaces_[surfi])
        {
            surfaceWriter& outWriter = writers_[surfi];

            if (outWriter.needsUpdate())
            {
                outWriter.setSurface(s);
            }

            outWriter.open(outputPath_/s.name());

            outWriter.beginTime(obr_.time());

            // Face-id field, but not if the writer manages that itself
            if (!s.isPointData() && s.hasFaceIds() && !outWriter.usesFaceIds())
            {
                // This is still quite horrible.

                Field<label> ids(s.faceIds());

                if
                (
                    returnReduce
                    (
                        !ListOps::found(ids, lessOp1<label>(0)),
                        andOp<bool>()
                    )
                )
                {
                    // From 0-based to 1-based, provided there are no
                    // negative ids (eg, encoded solid/side)

                    ids += label(1);
                }

                writeSurface(outWriter, ids, "Ids");
            }
        }
    }

    // Sample fields

    const IOobjectList objects(obr_, obr_.time().timeName());

    performAction<volScalarField>(objects, request);
    performAction<volVectorField>(objects, request);
    performAction<volSphericalTensorField>(objects, request);
    performAction<volSymmTensorField>(objects, request);
    performAction<volTensorField>(objects, request);

    // Only bother with surface fields if a sampler supports them
    if
    (
        testAny
        (
            surfaces(),
            [] (const sampledSurface& s) { return s.withSurfaceFields(); }
        )
    )
    {
        performAction<surfaceScalarField>(objects, request);
        performAction<surfaceVectorField>(objects, request);
        performAction<surfaceSphericalTensorField>(objects, request);
        performAction<surfaceSymmTensorField>(objects, request);
        performAction<surfaceTensorField>(objects, request);
    }


    // Finish this time step
    forAll(writers_, surfi)
    {
        if (((request & actions_[surfi]) & ACTION_WRITE) && nFaces_[surfi])
        {
            // Write geometry if no fields were written so that we still
            // can have something to look at

            if (!writers_[surfi].wroteData())
            {
                writers_[surfi].write();
            }

            writers_[surfi].endTime();
        }
    }

    return true;
}


bool Foam::sampledSurfaces::execute()
{
    if (onExecute_)
    {
        return performAction(ACTION_ALL & ~ACTION_WRITE);
    }

    return true;
}


bool Foam::sampledSurfaces::write()
{
    return performAction(ACTION_ALL);
}


void Foam::sampledSurfaces::updateMesh(const mapPolyMesh& mpm)
{
    if (&mpm.mesh() == &mesh_)
    {
        expire();
    }

    // pointMesh and interpolation will have been reset in mesh.update
}


void Foam::sampledSurfaces::movePoints(const polyMesh& mesh)
{
    if (&mesh == &mesh_)
    {
        expire();
    }
}


void Foam::sampledSurfaces::readUpdate(const polyMesh::readUpdateState state)
{
    if (state != polyMesh::UNCHANGED)
    {
        // May want to use force expiration here
        expire();
    }
}


bool Foam::sampledSurfaces::needsUpdate() const
{
    for (const sampledSurface& s : surfaces())
    {
        if (s.needsUpdate())
        {
            return true;
        }
    }

    return false;
}


bool Foam::sampledSurfaces::expire(const bool force)
{
    // Dimension as fraction of mesh bounding box
    const scalar mergeDim = mergeTol_ * mesh_.bounds().mag();

    label nChanged = 0;

    forAll(*this, surfi)
    {
        sampledSurface& s = (*this)[surfi];

        if (s.invariant() && !force)
        {
            // 'Invariant' - does not change when geometry does
            continue;
        }
        if (s.expire())
        {
            ++nChanged;
        }

        writers_[surfi].expire();
        writers_[surfi].mergeDim(mergeDim);
        nFaces_[surfi] = 0;
    }

    // True if any surfaces just expired
    return nChanged;
}


bool Foam::sampledSurfaces::update()
{
    if (!needsUpdate())
    {
        return false;
    }

    label nUpdated = 0;

    forAll(*this, surfi)
    {
        sampledSurface& s = (*this)[surfi];

        if (s.update())
        {
            ++nUpdated;
            writers_[surfi].expire();
        }

        nFaces_[surfi] = returnReduce(s.faces().size(), sumOp<label>());
    }

    return nUpdated;
}


Foam::scalar Foam::sampledSurfaces::mergeTol()
{
    return mergeTol_;
}


Foam::scalar Foam::sampledSurfaces::mergeTol(const scalar tol)
{
    const scalar prev(mergeTol_);
    mergeTol_ = tol;
    return prev;
}


// ************************************************************************* //
