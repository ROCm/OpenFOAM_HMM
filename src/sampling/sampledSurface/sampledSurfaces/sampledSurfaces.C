/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2018 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2017 OpenFOAM Foundation
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
#include "volFields.H"
#include "dictionary.H"
#include "Time.H"
#include "IOmanip.H"
#include "interpolationCell.H"
#include "volPointInterpolation.H"
#include "PatchTools.H"
#include "mapPolyMesh.H"
#include "sampledTriSurfaceMesh.H"
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

bool Foam::sampledSurfaces::verbose_ = false;
Foam::scalar Foam::sampledSurfaces::mergeTol_ = 1e-10;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sampledSurfaces::writeGeometry() const
{
    // Write to time directory under outputPath_
    // Skip surfaces without faces (eg, a failed cut-plane)

    const fileName outputDir = outputPath_/time_.timeName();

    forAll(*this, surfi)
    {
        const sampledSurface& s = operator[](surfi);

        if (Pstream::parRun())
        {
            if (Pstream::master() && mergedList_[surfi].size())
            {
                formatter_->write
                (
                    outputDir,
                    s.name(),
                    mergedList_[surfi]
                );
            }
        }
        else if (s.faces().size())
        {
            formatter_->write(outputDir, s.name(), s);
        }
    }
}


void Foam::sampledSurfaces::writeOriginalIds()
{
    const word fieldName = "Ids";
    const fileName outputDir = outputPath_/time_.timeName();

    forAll(*this, surfi)
    {
        const sampledSurface& s = operator[](surfi);

        if (s.hasFaceIds())
        {
            const labelList& idLst = s.originalIds();

            // Transcribe from label to scalar
            Field<scalar> ids(idLst.size());
            forAll(idLst, i)
            {
                ids[i] = idLst[i];
            }

            writeSurface(ids, surfi, fieldName, outputDir);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledSurfaces::sampledSurfaces
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObjects::regionFunctionObject(name, runTime, dict),
    PtrList<sampledSurface>(),
    mesh_(refCast<const fvMesh>(obr_)),
    loadFromFiles_(false),
    outputPath_
    (
        time_.globalPath()/functionObject::outputPrefix/name
    ),
    fieldSelection_(),
    sampleFaceScheme_(),
    sampleNodeScheme_(),
    mergedList_(),
    changedGeom_(),
    formatter_(nullptr)
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
    functionObjects::regionFunctionObject(name, obr, dict),
    PtrList<sampledSurface>(),
    mesh_(refCast<const fvMesh>(obr)),
    loadFromFiles_(loadFromFiles),
    outputPath_
    (
        time_.globalPath()/functionObject::outputPrefix/name
    ),
    fieldSelection_(),
    sampleFaceScheme_(),
    sampleNodeScheme_(),
    mergedList_(),
    changedGeom_(),
    formatter_(nullptr)
{
    outputPath_.clean();  // Remove unneeded ".."

    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::sampledSurfaces::verbose(const bool verbosity)
{
    verbose_ = verbosity;
}


bool Foam::sampledSurfaces::execute()
{
    return true;
}


bool Foam::sampledSurfaces::write()
{
    if (empty())
    {
        return true;
    }


    // Finalize surfaces, merge points etc.
    update();

    const label nFields = classifyFields();

    // Write geometry first if required,
    // or when no fields would otherwise be written
    if (formatter_->separateGeometry() || !nFields)
    {
        writeGeometry();
        changedGeom_ = false;
    }

    const IOobjectList objects(obr_, obr_.time().timeName());

    sampleAndWrite<volScalarField>(objects);
    sampleAndWrite<volVectorField>(objects);
    sampleAndWrite<volSphericalTensorField>(objects);
    sampleAndWrite<volSymmTensorField>(objects);
    sampleAndWrite<volTensorField>(objects);

    sampleAndWrite<surfaceScalarField>(objects);
    sampleAndWrite<surfaceVectorField>(objects);
    sampleAndWrite<surfaceSphericalTensorField>(objects);
    sampleAndWrite<surfaceSymmTensorField>(objects);
    sampleAndWrite<surfaceTensorField>(objects);

    return true;
}


bool Foam::sampledSurfaces::read(const dictionary& dict)
{
    PtrList<sampledSurface>::clear();
    mergedList_.clear();
    changedGeom_.clear();
    fieldSelection_.clear();

    sampleFaceScheme_ =
        dict.lookupOrDefault<word>("sampleScheme", "cell");

    sampleNodeScheme_ =
        dict.lookupOrDefault<word>("interpolationScheme", "cellPoint");

    const entry* eptr = dict.findEntry("surfaces");

    if (eptr && eptr->isDict())
    {
        PtrList<sampledSurface> surfs(eptr->dict().size());

        label surfi = 0;

        for (const entry& dEntry : eptr->dict())
        {
            if (dEntry.isDict())
            {
                autoPtr<sampledSurface> surf =
                    sampledSurface::New
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
        PtrList<sampledSurface> surfs
        (
            eptr->stream(),
            sampledSurface::iNew(mesh_)
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


    const auto& surfs = surfaces();
    if (surfs.size())
    {
        dict.readEntry("fields", fieldSelection_);
        fieldSelection_.uniq();

        // The surface writer and format options
        const word writerType(dict.get<word>("surfaceFormat"));
        const dictionary& formatOptions = dict.subOrEmptyDict("formatOptions");

        // Define the surface formatter
        // Optionally defined extra controls for the output formats
        formatter_ = surfaceWriter::New
        (
            writerType,
            dict.subOrEmptyDict("formatOptions").subOrEmptyDict(writerType)
        );

        if (Pstream::parRun())
        {
            mergedList_.resize(size());
        }

        // Ensure all surfaces and merge information are expired
        expire();

        label surfi = 0;
        for (const sampledSurface& s : surfs)
        {
            if (!surfi)
            {
                Info<< "Sampled surface:" << nl;
            }
            Info<< "    " << s.name() << " -> " << writerType << nl;

            ++surfi;
        }
        Info<< nl;
    }

    if (debug && Pstream::master())
    {
        Pout<< "sample fields:" << fieldSelection_ << nl
            << "sample surfaces:" << nl << "(" << nl;

        for (const sampledSurface& s : surfaces())
        {
            Pout<< "  " << s << nl;
        }

        Pout<< ")" << endl;
    }

    // New geometry
    changedGeom_.resize(size());
    changedGeom_ = true;

    return true;
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
        expire();
    }
}


bool Foam::sampledSurfaces::needsUpdate() const
{
    forAll(*this, surfi)
    {
        if (operator[](surfi).needsUpdate())
        {
            return true;
        }
    }

    return false;
}


bool Foam::sampledSurfaces::expire()
{
    bool justExpired = false;

    forAll(*this, surfi)
    {
        if (operator[](surfi).expire())
        {
            justExpired = true;
        }

        // Clear merge information
        if (Pstream::parRun())
        {
            mergedList_[surfi].clear();
        }
    }

    changedGeom_ = true;

    // true if any surfaces just expired
    return justExpired;
}


bool Foam::sampledSurfaces::update()
{
    if (!needsUpdate())
    {
        return false;
    }

    bool updated = false;

    // Serial: quick and easy, no merging required
    if (!Pstream::parRun())
    {
        forAll(*this, surfi)
        {
            sampledSurface& s = operator[](surfi);

            if (s.update())
            {
                updated = true;
                changedGeom_[surfi] = true;
            }
        }

        return updated;
    }


    // Dimension as fraction of mesh bounding box
    const scalar mergeDim = mergeTol_*mesh_.bounds().mag();

    if (Pstream::master() && debug)
    {
        Pout<< nl << "Merging all points within "
            << mergeDim << " metre" << endl;
    }

    forAll(*this, surfi)
    {
        sampledSurface& s = operator[](surfi);

        if (s.update())
        {
            updated = true;
            changedGeom_[surfi] = true;
            mergedList_[surfi].merge(s, mergeDim);
        }
    }

    return updated;
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
