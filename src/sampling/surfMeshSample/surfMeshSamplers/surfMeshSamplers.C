/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2019 OpenCFD Ltd.
     \\/     M anipulation  |
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

#include "surfMeshSamplers.H"
#include "volFields.H"
#include "dictionary.H"
#include "Time.H"
#include "IOmanip.H"
#include "volPointInterpolation.H"
#include "PatchTools.H"
#include "mapPolyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(surfMeshSamplers, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        surfMeshSamplers,
        dictionary
    );
} // End namespace Foam


bool Foam::surfMeshSamplers::verbose_ = false;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// //- Temporary object registry for passing around values
// const objectRegistry& tmpRegistry() const;
//
// //- Temporary object registry for passing around values
// const objectRegistry& tmpRegistry(const word& subName) const;

// const Foam::objectRegistry&
// Foam::surfMeshSamplers::tmpRegistry() const
// {
//     // Sub-registry for sampling, choose name for fewer collisions
//     return mesh_.thisDb().subRegistry
//     (
//         "$tmp$" + type() + "$" + name(),
//         true,
//         false
//     );
// }
//
//
// const Foam::objectRegistry&
// Foam::surfMeshSamplers::tmpRegistry(const word& subName) const
// {
//     return tmpRegistry().subRegistry(subName, true, false);
// }


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfMeshSamplers::surfMeshSamplers
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObjects::fvMeshFunctionObject(name, runTime, dict),
    PtrList<surfMeshSample>(),
    fieldSelection_(),
    sampleFaceScheme_()
{
    read(dict);
}


Foam::surfMeshSamplers::surfMeshSamplers
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    functionObjects::fvMeshFunctionObject(name, obr, dict),
    PtrList<surfMeshSample>(),
    fieldSelection_(),
    sampleFaceScheme_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::surfMeshSamplers::verbose(const bool verbosity)
{
    verbose_ = verbosity;
}


bool Foam::surfMeshSamplers::execute()
{
    if (empty())
    {
        return true;
    }

    // The acceptable fields
    wordHashSet acceptable;
    acceptable.insert(acceptType<scalar>());
    acceptable.insert(acceptType<vector>());
    acceptable.insert(acceptType<sphericalTensor>());
    acceptable.insert(acceptType<symmTensor>());
    acceptable.insert(acceptType<tensor>());

    const wordList fields = acceptable.sortedToc();
    if (!fields.empty())
    {
        for (surfMeshSample& s : surfaces())
        {
            // Potentially monitor the update for writing geometry?
            if (s.needsUpdate())
            {
                s.update();
            }

            s.sample(fields, sampleFaceScheme_);
        }
    }

    clearObjects(removeFieldsOnExecute_);

    return true;
}


bool Foam::surfMeshSamplers::write()
{
    // Write sampled fields (on surface)
    //
    // Doesn't bother checking which fields have been generated here
    // or elsewhere

    for (const surfMeshSample& s : surfaces())
    {
        s.write(fieldSelection_);
    }

    clearObjects(removeFieldsOnWrite_);

    return true;
}


bool Foam::surfMeshSamplers::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    PtrList<surfMeshSample>::clear();
    fieldSelection_.clear();
    removeFieldsOnExecute_.clear();
    removeFieldsOnWrite_.clear();

    dict.readIfPresent("removeFieldsOnExecute", removeFieldsOnExecute_);
    dict.readIfPresent("removeFieldsOnWrite", removeFieldsOnWrite_);

    const bool createOnRead = dict.lookupOrDefault("createOnRead", false);

    sampleFaceScheme_ =
        dict.lookupOrDefault<word>("sampleScheme", "cell");

    const entry* eptr = dict.findEntry("surfaces");

    if (eptr && eptr->isDict())
    {
        PtrList<surfMeshSample> surfs(eptr->dict().size());

        label surfi = 0;

        for (const entry& dEntry : eptr->dict())
        {
            if (dEntry.isDict())
            {
                autoPtr<surfMeshSample> surf =
                    surfMeshSample::New
                    (
                        dEntry.keyword(),
                        mesh_,
                        dEntry.dict()
                    );

                if (surf.valid() && surf->enabled())
                {
                    surfs.set(surfi, surf);
                    ++surfi;
                }
            }
        }

        surfs.resize(surfi);
        surfaces().transfer(surfs);
    }
    else if (eptr)
    {
        PtrList<surfMeshSample> surfs
        (
            eptr->stream(),
            surfMeshSample::iNew(mesh_)
        );

        forAll(surfs, surfi)
        {
            if (!surfs[surfi].enabled())
            {
                surfs.set(surfi, nullptr);
            }
        }

        surfs.resize(surfs.squeezeNull());

        surfaces().transfer(surfs);
    }


    auto& surfs = surfaces();
    if (surfs.size())
    {
        dict.readEntry("fields", fieldSelection_);
        fieldSelection_.uniq();

        Info<< type() << " fields:  " << flatOutput(fieldSelection_) << nl;
        Info<< nl;

        // Need to initialize corresponding surfMesh for others in the chain.
        // This can simply be a zero-sized placeholder, or the real surface with
        // faces.

        label surfi = 0;
        for (surfMeshSample& s : surfs)
        {
            if (!surfi)
            {
                Info<< "Sampled surface:" << nl;
            }
            Info<< "    " << s.name() << nl;

            if (createOnRead)
            {
                s.update();
            }
            else
            {
                s.create();
            }

            ++surfi;
        }

        Info<< nl;
    }

    if (removeFieldsOnExecute_.size())
    {
        Info<< type() << " Remove fields on-execute : "
            << removeFieldsOnExecute_ << nl;
    }
    if (removeFieldsOnWrite_.size())
    {
        Info<< type() << " Remove fields on-write   :"
            << removeFieldsOnWrite_ << nl;
    }

    // Ensure all surfaces are expired (unsampled)
    expire();

    return true;
}


void Foam::surfMeshSamplers::updateMesh(const mapPolyMesh& mpm)
{
    if (&mpm.mesh() == &mesh_)
    {
        expire();
    }

    // pointMesh and interpolation will have been reset in mesh.update
}


void Foam::surfMeshSamplers::movePoints(const polyMesh& m)
{
    if (&m == &mesh_)
    {
        expire();
    }
}


void Foam::surfMeshSamplers::readUpdate(const polyMesh::readUpdateState state)
{
    if (state != polyMesh::UNCHANGED)
    {
        expire();
    }
}


bool Foam::surfMeshSamplers::needsUpdate() const
{
    for (const surfMeshSample& s : surfaces())
    {
        if (s.needsUpdate())
        {
            return true;
        }
    }

    return false;
}


bool Foam::surfMeshSamplers::expire()
{
    label nChanged = 0;

    for (surfMeshSample& s : surfaces())
    {
        if (s.expire())
        {
            ++nChanged;
        }
    }

    // True if any surfaces just expired
    return nChanged;
}


bool Foam::surfMeshSamplers::update()
{
    if (!needsUpdate())
    {
        return false;
    }

    label nUpdated = 0;

    for (surfMeshSample& s : surfaces())
    {
        if (s.update())
        {
            ++nUpdated;
        }
    }

    return nUpdated;
}


// ************************************************************************* //
