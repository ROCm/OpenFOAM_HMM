/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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
#include "wordReListMatcher.H"
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
}


bool Foam::surfMeshSamplers::verbose_ = false;

void Foam::surfMeshSamplers::checkOutNames
(
    const objectRegistry& registry,
    const UList<word>& names
)
{
    objectRegistry& reg = const_cast<objectRegistry&>(registry);

    forAll(names, namei)
    {
        objectRegistry::iterator iter = reg.find(names[namei]);
        if (iter != reg.end())
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfMeshSamplers::surfMeshSamplers
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObjects::regionFunctionObject(name, runTime, dict),
    PtrList<surfMeshSampler>(),
    mesh_(refCast<const fvMesh>(obr_)),
    fieldSelection_()
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
    PtrList<surfMeshSampler>(),
    mesh_(refCast<const fvMesh>(obr)),
    fieldSelection_(),
    derivedNames_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfMeshSamplers::~surfMeshSamplers()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::surfMeshSamplers::verbose(const bool verbosity)
{
    verbose_ = verbosity;
}


bool Foam::surfMeshSamplers::execute()
{
    if (size())
    {
        const objectRegistry& db = mesh_.thisDb();

        // Manage derived names
        DynamicList<word> added(derivedNames_.size());
        DynamicList<word> cleanup(derivedNames_.size());

        forAll(derivedNames_, namei)
        {
            const word& derivedName = derivedNames_[namei];

            if (derivedName == "rhoU")
            {
                added.append(derivedName);

                if (!db.foundObject<volVectorField>(derivedName))
                {
                    cleanup.append(derivedName);

                    db.store
                    (
                        new volVectorField
                        (
                            derivedName,
                            // rhoU = rho * U
                            (
                                mesh_.lookupObject<volScalarField>("rho")
                              * mesh_.lookupObject<volVectorField>("U")
                            )
                        )
                    );
                }
            }
            else if (derivedName == "pTotal")
            {
                added.append(derivedName);

                if (!db.foundObject<volScalarField>(derivedName))
                {
                    cleanup.append(derivedName);

                    db.store
                    (
                        new volScalarField
                        (
                            derivedName,
                            // pTotal = p + U^2 / 2
                            (
                                mesh_.lookupObject<volScalarField>("p")
                              + 0.5
                              * mesh_.lookupObject<volScalarField>("rho")
                              * magSqr(mesh_.lookupObject<volVectorField>("U"))
                            )
                        )
                    );
                }
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
            forAll(*this, surfI)
            {
                surfMeshSampler& s = operator[](surfI);

                // Potentially monitor the update for writing geometry?
                if (s.needsUpdate())
                {
                    s.update();
                }

                s.sample(fields);
            }
        }

        checkOutNames(db, cleanup);
    }

    return true;
}


bool Foam::surfMeshSamplers::write()
{
    // Write sampled fields (on surface)
    //
    // Doesn't bother checking which fields have been generated here
    // or elsewhere

    // This could be more efficient
    wordReList select(fieldSelection_.size() + derivedNames_.size());

    label nElem = 0;
    forAll(fieldSelection_, i)
    {
        select[nElem++] = fieldSelection_[i];
    }
    forAll(derivedNames_, i)
    {
        select[nElem++] = derivedNames_[i];
    }

    // avoid duplicate entries
    select = wordReListMatcher::uniq(select);

    forAll(*this, surfI)
    {
        const surfMeshSampler& s = operator[](surfI);
        s.write(select);
    }

    return true;
}


bool Foam::surfMeshSamplers::read(const dictionary& dict)
{
    derivedNames_.clear();

    const bool createOnRead =
        dict.lookupOrDefault<Switch>("createOnRead", false);

    if (dict.found("surfaces"))
    {
        fieldSelection_ = wordReListMatcher::uniq
        (
            wordReList(dict.lookup("fields"))
        );
        Info<< type() << " fields: " << fieldSelection_ << nl;

        if (dict.readIfPresent("derived", derivedNames_))
        {
            Info<< type() << " derived: " << derivedNames_ << nl;
        }

        PtrList<surfMeshSampler> newList
        (
            dict.lookup("surfaces"),
            surfMeshSampler::iNew(mesh_)
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
            forAll(*this, surfI)
            {
                surfMeshSampler& s = operator[](surfI);

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
    forAll(*this, surfI)
    {
        if (operator[](surfI).needsUpdate())
        {
            return true;
        }
    }

    return false;
}


bool Foam::surfMeshSamplers::expire()
{
    bool justExpired = false;

    forAll(*this, surfI)
    {
        if (operator[](surfI).expire())
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
    forAll(*this, surfI)
    {
        if (operator[](surfI).update())
        {
            updated = true;
        }
    }

    return updated;
}


// ************************************************************************* //
