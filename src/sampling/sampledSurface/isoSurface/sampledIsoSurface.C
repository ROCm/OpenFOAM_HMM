/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2019 OpenCFD Ltd.
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

#include "sampledIsoSurface.H"
#include "dictionary.H"
#include "volFields.H"
#include "volPointInterpolation.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledIsoSurface, 0);
    addNamedToRunTimeSelectionTable
    (
        sampledSurface,
        sampledIsoSurface,
        word,
        isoSurface
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sampledIsoSurface::getIsoFields() const
{
    const fvMesh& fvm = static_cast<const fvMesh&>(mesh());

    // Get volField
    // ~~~~~~~~~~~~

    volFieldPtr_ = fvm.getObjectPtr<volScalarField>(isoField_);

    if (volFieldPtr_)
    {
        DebugInFunction
            << "Lookup volField " << isoField_ << endl;

        storedVolFieldPtr_.clear();
    }
    else
    {
        // Bit of a hack. Read field and store.

        DebugInFunction
            << "Checking " << isoField_
            << " for same time " << fvm.time().timeName() << endl;

        if
        (
            storedVolFieldPtr_.empty()
         || (fvm.time().timeName() != storedVolFieldPtr_().instance())
        )
        {
            DebugInFunction
                << "Reading volField " << isoField_
                << " from time " << fvm.time().timeName() << endl;

            IOobject vfHeader
            (
                isoField_,
                fvm.time().timeName(),
                fvm,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            );

            if (vfHeader.typeHeaderOk<volScalarField>(true))
            {
                storedVolFieldPtr_.reset
                (
                    new volScalarField
                    (
                        vfHeader,
                        fvm
                    )
                );
                volFieldPtr_ = storedVolFieldPtr_.get(); // get(), not release()
            }
            else
            {
                FatalErrorInFunction
                    << "Cannot find isosurface field " << isoField_
                    << " in database or directory " << vfHeader.path()
                    << exit(FatalError);
            }
        }
    }


    // Get pointField
    // ~~~~~~~~~~~~~~

    // In case of multiple iso values we don't want to calculate multiple e.g.
    // "volPointInterpolate(p)" so register it and re-use it. This is the
    // same as the 'cache' functionality from volPointInterpolate but
    // unfortunately that one does not guarantee that the field pointer
    // remain: e.g. some other functionObject might delete the cached version.
    // (volPointInterpolation::interpolate with cache=false deletes any
    //  registered one or if mesh.changing())

    if (subMeshPtr_.empty())
    {
        const word pointFldName =
            "volPointInterpolate_"
          + type()
          + "("
          + isoField_
          + ')';


        pointFieldPtr_ = fvm.getObjectPtr<pointScalarField>(pointFldName);

        if (pointFieldPtr_)
        {
            DebugInFunction
                << "lookup pointField " << pointFldName << endl;

            if (!pointFieldPtr_->upToDate(*volFieldPtr_))
            {
                DebugInFunction
                    << "updating pointField " << pointFldName << endl;

                // Update the interpolated value
                volPointInterpolation::New(fvm).interpolate
                (
                    *volFieldPtr_,
                    const_cast<pointScalarField&>(*pointFieldPtr_)
                );
            }
        }
        else
        {
            // Not in registry. Interpolate.

            DebugInFunction
                << "creating pointField " << pointFldName << endl;

            // Interpolate without cache. Note that we're registering it
            // below so next time round it goes into the condition
            // above.
            pointFieldPtr_ =
                volPointInterpolation::New(fvm).interpolate
                (
                    *volFieldPtr_,
                    pointFldName,
                    false
                ).ptr();

            const_cast<pointScalarField*>(pointFieldPtr_)->store();
        }


        // If averaging redo the volField.
        // Can only be done now since needs the point field.
        if (average_)
        {
            storedVolFieldPtr_.reset
            (
                pointAverage(*pointFieldPtr_).ptr()
            );
            volFieldPtr_ = storedVolFieldPtr_.get(); // get(), not release()
        }


        DebugInFunction
            << "volField " << volFieldPtr_->name()
            << " min:" << min(*volFieldPtr_).value()
            << " max:" << max(*volFieldPtr_).value() << nl
            << "pointField " << pointFieldPtr_->name()
            << " min:" << gMin(pointFieldPtr_->primitiveField())
            << " max:" << gMax(pointFieldPtr_->primitiveField()) << endl;
    }
    else
    {
        // Get subMesh variants
        const fvMesh& subFvm = subMeshPtr_().subMesh();

        // Either lookup on the submesh or subset the whole-mesh volField

        volSubFieldPtr_ = subFvm.getObjectPtr<volScalarField>(isoField_);

        if (volSubFieldPtr_)
        {
            DebugInFunction
                << "Sub-mesh lookup volField " << isoField_ << endl;

            storedVolSubFieldPtr_.clear();
        }
        else
        {
            DebugInFunction
                << "Sub-setting volField " << isoField_ << endl;

            storedVolSubFieldPtr_.reset
            (
                subMeshPtr_().interpolate
                (
                    *volFieldPtr_
                ).ptr()
            );
            storedVolSubFieldPtr_->checkOut();
            volSubFieldPtr_ = storedVolSubFieldPtr_.get(); // not release()
        }


        // The point field on subMesh

        const word pointFldName =
            "volPointInterpolate_"
          + type()
          + "("
          + volSubFieldPtr_->name()
          + ')';


        pointFieldPtr_ = subFvm.getObjectPtr<pointScalarField>(pointFldName);

        if (pointFieldPtr_)
        {
            DebugInFunction
                << "Sub-mesh lookup pointField " << pointFldName << endl;

            if (!pointFieldPtr_->upToDate(*volSubFieldPtr_))
            {
                DebugInFunction
                    << "Updating submesh pointField " << pointFldName << endl;

                // Update the interpolated value
                volPointInterpolation::New(subFvm).interpolate
                (
                    *volSubFieldPtr_,
                    const_cast<pointScalarField&>(*pointFieldPtr_)
                );
            }
        }
        else
        {
            DebugInFunction
                << "Interpolating submesh volField "
                << volSubFieldPtr_->name()
                << " to get submesh pointField " << pointFldName << endl;

            pointSubFieldPtr_ =
                volPointInterpolation::New
                (
                    subFvm
                ).interpolate(*volSubFieldPtr_).ptr();

            const_cast<pointScalarField*>(pointSubFieldPtr_)->store();
        }


        // If averaging redo the volField. Can only be done now since needs the
        // point field.
        if (average_)
        {
            storedVolSubFieldPtr_.reset
            (
                pointAverage(*pointSubFieldPtr_).ptr()
            );
            volSubFieldPtr_ = storedVolSubFieldPtr_.get(); // not release()
        }


        DebugInFunction
            << "volSubField "
            << volSubFieldPtr_->name()
            << " min:" << min(*volSubFieldPtr_).value()
            << " max:" << max(*volSubFieldPtr_).value() << nl
            << "pointSubField "
            << pointSubFieldPtr_->name()
            << " min:" << gMin(pointSubFieldPtr_->primitiveField())
            << " max:" << gMax(pointSubFieldPtr_->primitiveField()) << endl;
    }
}


bool Foam::sampledIsoSurface::updateGeometry() const
{
    const fvMesh& fvm = static_cast<const fvMesh&>(mesh());

    // No update needed
    if (fvm.time().timeIndex() == prevTimeIndex_)
    {
        return false;
    }

    // Get sub-mesh if any
    if
    (
        (-1 != mesh().cellZones().findIndex(zoneNames_))
     && subMeshPtr_.empty()
    )
    {
        const polyBoundaryMesh& patches = mesh().boundaryMesh();

        // Patch to put exposed internal faces into
        const label exposedPatchi = patches.findPatchID(exposedPatchName_);

        DebugInfo
            << "Allocating subset of size "
            << mesh().cellZones().selection(zoneNames_).count()
            << " with exposed faces into patch "
            << patches[exposedPatchi].name() << endl;

        subMeshPtr_.reset
        (
            new fvMeshSubset
            (
                fvm,
                mesh().cellZones().selection(zoneNames_),
                exposedPatchi
            )
        );
    }


    prevTimeIndex_ = fvm.time().timeIndex();
    getIsoFields();

    // Clear any stored topo
    surfPtr_.clear();

    // Clear derived data
    clearGeom();

    if (subMeshPtr_.valid())
    {
        const volScalarField& vfld = *volSubFieldPtr_;

        surfPtr_.reset
        (
            new isoSurface
            (
                vfld,
                *pointSubFieldPtr_,
                isoVal_,
                filter_,
                bounds_,
                mergeTol_
            )
        );
    }
    else
    {
        const volScalarField& vfld = *volFieldPtr_;

        surfPtr_.reset
        (
            new isoSurface
            (
                vfld,
                *pointFieldPtr_,
                isoVal_,
                filter_,
                bounds_,
                mergeTol_
            )
        );
    }


    if (debug)
    {
        Pout<< "sampledIsoSurface::updateGeometry() : constructed iso:"
            << nl
            << "    filter         : " << Switch(bool(filter_)) << nl
            << "    average        : " << Switch(average_) << nl
            << "    isoField       : " << isoField_ << nl
            << "    isoValue       : " << isoVal_ << nl;
        if (subMeshPtr_.valid())
        {
            Pout<< "    zone size      : " << subMeshPtr_().subMesh().nCells()
                << nl;
        }
        Pout<< "    points         : " << points().size() << nl
            << "    faces          : " << surface().size() << nl
            << "    cut cells      : " << surface().meshCells().size()
            << endl;
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledIsoSurface::sampledIsoSurface
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledSurface(name, mesh, dict),
    isoField_(dict.get<word>("isoField")),
    isoVal_(dict.get<scalar>("isoValue")),
    mergeTol_(dict.getOrDefault<scalar>("mergeTol", 1e-6)),
    filter_
    (
        isoSurfaceBase::getFilterType
        (
            dict,
            isoSurfaceBase::filterType::DIAGCELL
        )
    ),
    average_(dict.getOrDefault("average", false)),
    bounds_(dict.getOrDefault("bounds", boundBox::invertedBox)),
    zoneNames_(),
    exposedPatchName_(),
    surfPtr_(nullptr),
    prevTimeIndex_(-1),
    storedVolFieldPtr_(nullptr),
    volFieldPtr_(nullptr),
    pointFieldPtr_(nullptr)
{
    if (!sampledSurface::interpolate())
    {
        FatalIOErrorInFunction(dict)
            << "Non-interpolated iso surface not supported since triangles"
            << " span across cells." << exit(FatalIOError);
    }


    if (!dict.readIfPresent("zones", zoneNames_) && dict.found("zone"))
    {
        zoneNames_.resize(1);
        dict.readEntry("zone", zoneNames_.first());
    }

    if (-1 != mesh.cellZones().findIndex(zoneNames_))
    {
        dict.readEntry("exposedPatchName", exposedPatchName_);

        if (mesh.boundaryMesh().findPatchID(exposedPatchName_) == -1)
        {
            FatalIOErrorInFunction(dict)
                << "Cannot find patch " << exposedPatchName_
                << " in which to put exposed faces." << endl
                << "Valid patches are " << mesh.boundaryMesh().names()
                << exit(FatalIOError);
        }

        DebugInfo
            << "Restricting to cellZone(s) " << flatOutput(zoneNames_)
            << " with exposed internal faces into patch "
            << exposedPatchName_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sampledIsoSurface::~sampledIsoSurface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sampledIsoSurface::needsUpdate() const
{
    const fvMesh& fvm = static_cast<const fvMesh&>(mesh());

    return fvm.time().timeIndex() != prevTimeIndex_;
}


bool Foam::sampledIsoSurface::expire()
{
    surfPtr_.clear();
    subMeshPtr_.clear();

    // Clear derived data
    clearGeom();

    // Already marked as expired
    if (prevTimeIndex_ == -1)
    {
        return false;
    }

    // Force update
    prevTimeIndex_ = -1;
    return true;
}


bool Foam::sampledIsoSurface::update()
{
    return updateGeometry();
}


Foam::tmp<Foam::scalarField> Foam::sampledIsoSurface::sample
(
    const interpolation<scalar>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::vectorField> Foam::sampledIsoSurface::sample
(
    const interpolation<vector>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::sphericalTensorField> Foam::sampledIsoSurface::sample
(
    const interpolation<sphericalTensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledIsoSurface::sample
(
    const interpolation<symmTensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::tensorField> Foam::sampledIsoSurface::sample
(
    const interpolation<tensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::scalarField> Foam::sampledIsoSurface::interpolate
(
    const interpolation<scalar>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::vectorField> Foam::sampledIsoSurface::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}

Foam::tmp<Foam::sphericalTensorField> Foam::sampledIsoSurface::interpolate
(
    const interpolation<sphericalTensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::symmTensorField> Foam::sampledIsoSurface::interpolate
(
    const interpolation<symmTensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::tensorField> Foam::sampledIsoSurface::interpolate
(
    const interpolation<tensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


void Foam::sampledIsoSurface::print(Ostream& os) const
{
    os  << "sampledIsoSurface: " << name() << " :"
        << "  field   :" << isoField_
        << "  value   :" << isoVal_;
}


// ************************************************************************* //
