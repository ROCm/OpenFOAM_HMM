/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2018 OpenCFD Ltd.
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

void Foam::surfMeshSamplers::checkOutNames
(
    const objectRegistry& registry,
    const UList<word>& names
)
{
    objectRegistry& reg = const_cast<objectRegistry&>(registry);

    for (const word& fldName : names)
    {
        objectRegistry::iterator iter = reg.find(fldName);
        if (iter.found())
        {
            registry.checkOut(*iter());
        }
    }
}


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

bool Foam::surfMeshSamplers::add_rhoU(const word& derivedName)
{
    const objectRegistry& db = mesh_.thisDb();

    const bool existed = db.foundObject<volVectorField>(derivedName);

    if (existed)
    {
        return false; // Volume field already existed - not added.
    }


    // rhoU = rho * U

    const auto* rhoPtr = mesh_.findObject<volScalarField>("rho");
    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");

    tmp<volVectorField> tresult;

    if (rhoPtr)
    {
        const auto& rho = *rhoPtr;

        tresult = tmp<volVectorField>::New
        (
            derivedName, (rho * U)
        );
    }
    else
    {
        const dimensionedScalar rho("rho", dimDensity, rhoRef_);

        tresult = tmp<volVectorField>::New
        (
            derivedName, (rho * U)
        );
    }

    db.store(tresult.ptr());

    return !existed;
}


bool Foam::surfMeshSamplers::add_pTotal(const word& derivedName)
{
    const objectRegistry& db = mesh_.thisDb();

    const bool existed = db.foundObject<volVectorField>(derivedName);

    if (existed)
    {
        return false; // Volume field already existed - not added.
    }

    // pTotal = p + rho * U^2 / 2

    const auto* rhoPtr = mesh_.findObject<volScalarField>("rho");
    const volScalarField& p = mesh_.lookupObject<volScalarField>("p");
    const volVectorField& U = mesh_.lookupObject<volVectorField>("U");

    tmp<volScalarField> tresult;

    if (rhoPtr)
    {
        const auto& rho = *rhoPtr;

        tresult = tmp<volScalarField>::New
        (
            derivedName, (p + 0.5 * rho * magSqr(U))
        );
    }
    else
    {
        const dimensionedScalar rho("rho", dimDensity, rhoRef_);

        tresult = tmp<volScalarField>::New
        (
            derivedName, (rho * (p + 0.5 * magSqr(U)))
        );
    }

    db.store(tresult.ptr());

    return !existed;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfMeshSamplers::surfMeshSamplers
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObjects::regionFunctionObject(name, runTime, dict),
    PtrList<surfMeshSample>(),
    mesh_(refCast<const fvMesh>(obr_)),
    fieldSelection_(),
    derivedNames_(),
    sampleScheme_(word::null),
    rhoRef_(1.0)
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
    functionObjects::regionFunctionObject(name, obr, dict),
    PtrList<surfMeshSample>(),
    mesh_(refCast<const fvMesh>(obr)),
    fieldSelection_(),
    derivedNames_(),
    sampleScheme_(word::null),
    rhoRef_(1.0)
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

    const objectRegistry& db = mesh_.thisDb();

    // Manage derived names
    DynamicList<word> added(derivedNames_.size());
    DynamicList<word> cleanup(derivedNames_.size());

    for (const word& derivedName : derivedNames_)
    {
        // This is a fairly ugly dispatch mechanism

        if (derivedName == "rhoU")
        {
            if (add_rhoU(derivedName))
            {
                cleanup.append(derivedName);
            }
            added.append(derivedName);
        }
        else if (derivedName == "pTotal")
        {
            if (add_pTotal(derivedName))
            {
                cleanup.append(derivedName);
            }
            added.append(derivedName);
        }
        else
        {
            WarningInFunction
                << "unknown derived name: " << derivedName << nl
                << "Use one of 'rhoU', 'pTotal'" << nl
                << endl;
        }
    }


    // The acceptable fields
    wordHashSet acceptable(added);
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

            s.sample(fields, sampleScheme_);
        }
    }

    checkOutNames(db, cleanup);

    return true;
}


bool Foam::surfMeshSamplers::write()
{
    // Write sampled fields (on surface)
    //
    // Doesn't bother checking which fields have been generated here
    // or elsewhere

    // This could be more efficient
    wordRes select(fieldSelection_.size() + derivedNames_.size());

    label nElem = 0;
    for (const wordRe& item : fieldSelection_)
    {
        select[nElem++] = item;
    }
    for (const auto& derivedName : derivedNames_)
    {
        select[nElem++] = derivedName;
    }

    // Avoid duplicate entries
    select.uniq();

    for (const surfMeshSample& s : surfaces())
    {
        s.write(select);
    }

    return true;
}


bool Foam::surfMeshSamplers::read(const dictionary& dict)
{
    fieldSelection_.clear();
    derivedNames_.clear();

    rhoRef_ = dict.lookupOrDefault<scalar>("rhoRef", 1);

    const bool createOnRead = dict.lookupOrDefault("createOnRead", false);

    if (dict.found("surfaces"))
    {
        sampleScheme_ = dict.lookupOrDefault<word>("sampleScheme", "cell");

        dict.readEntry("fields", fieldSelection_);
        fieldSelection_.uniq();

        Info<< type() << " fields:  " << flatOutput(fieldSelection_) << nl;

        if (dict.readIfPresent("derived", derivedNames_))
        {
            Info<< type() << " derived: " << flatOutput(derivedNames_) << nl;
        }
        Info<< nl;

        PtrList<surfMeshSample> newList
        (
            dict.lookup("surfaces"),
            surfMeshSample::iNew(mesh_)
        );
        transfer(newList);

        // Ensure all surfaces and merge information are expired
        expire();

        // Need to initialize corresponding surfMesh for others in the chain.
        // This can simply be a zero-sized placeholder, or the real surface with
        // faces.
        if (this->size())
        {
            Info<< "Reading surface description:" << nl;
            for (surfMeshSample& s : surfaces())
            {
                Info<< "    " << s.name() << nl;
                if (createOnRead)
                {
                    s.update();
                }
                else
                {
                    s.create();
                }
            }
            Info<< endl;
        }
    }

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
    bool justExpired = false;

    for (surfMeshSample& s : surfaces())
    {
        if (s.expire())
        {
            justExpired = true;
        }
    }

    // True if any surfaces just expired
    return justExpired;
}


bool Foam::surfMeshSamplers::update()
{
    if (!needsUpdate())
    {
        return false;
    }

    bool updated = false;
    for (surfMeshSample& s : surfaces())
    {
        if (s.update())
        {
            updated = true;
        }
    }

    return updated;
}


// ************************************************************************* //
