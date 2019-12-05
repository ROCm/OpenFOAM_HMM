/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "sampledCuttingPlane.H"
#include "dictionary.H"
#include "volFields.H"
#include "volPointInterpolation.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sampledCuttingPlane, 0);
    addNamedToRunTimeSelectionTable
    (
        sampledSurface,
        sampledCuttingPlane,
        word,
        cuttingPlane
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::sampledCuttingPlane::checkBoundsIntersection
(
    const plane& pln,
    const boundBox& meshBb
) const
{
    // Verify specified bounding box
    if (bounds_.valid())
    {
        // Bounding box does not overlap with (global) mesh!
        if (!bounds_.overlaps(meshBb))
        {
            WarningInFunction
                << nl
                << name() << " : "
                << "Bounds " << bounds_
                << " do not overlap the mesh bounding box " << meshBb
                << nl << endl;
        }

        // Plane does not intersect the bounding box
        if (!bounds_.intersects(pln))
        {
            WarningInFunction
                << nl
                << name() << " : "
                << "Plane "<< pln << " does not intersect the bounds "
                << bounds_
                << nl << endl;
        }
    }

    // Plane does not intersect the (global) mesh!
    if (!meshBb.intersects(pln))
    {
        WarningInFunction
            << nl
            << name() << " : "
            << "Plane "<< pln << " does not intersect the mesh bounds "
            << meshBb
            << nl << endl;
    }
}


void Foam::sampledCuttingPlane::createGeometry()
{
    if (debug)
    {
        Pout<< "sampledCuttingPlane::createGeometry :updating geometry."
            << endl;
    }

    // Clear any stored topologies
    isoSurfPtr_.clear();
    isoSurfCellPtr_.clear();
    isoSurfTopoPtr_.clear();
    pointDistance_.clear();
    cellDistancePtr_.clear();

    // Clear derived data
    clearGeom();

    const fvMesh& fvm = static_cast<const fvMesh&>(this->mesh());

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

        bitSet cellsToSelect = mesh().cellZones().selection(zoneNames_);

        DebugInfo
            << "Allocating subset of size "
            << cellsToSelect.count()
            << " with exposed faces into patch "
            << patches[exposedPatchi].name() << endl;


        // If we will use a fvMeshSubset so can apply bounds as well to make
        // the initial selection smaller.
        if (bounds_.valid() && cellsToSelect.any())
        {
            const auto& cellCentres = fvm.C();

            for (const label celli : cellsToSelect)
            {
                const point& cc = cellCentres[celli];

                if (!bounds_.contains(cc))
                {
                    cellsToSelect.unset(celli);
                }
            }

            DebugInfo
                << "Bounded subset of size "
                << cellsToSelect.count() << endl;
        }

        subMeshPtr_.reset(new fvMeshSubset(fvm, cellsToSelect, exposedPatchi));
    }


    // Select either the submesh or the underlying mesh
    const fvMesh& mesh =
    (
        subMeshPtr_.valid()
      ? subMeshPtr_().subMesh()
      : fvm
    );

    checkBoundsIntersection(plane_, mesh.bounds());


    // Distance to cell centres
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    cellDistancePtr_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "cellDistance",
                mesh.time().timeName(),
                mesh.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimensionedScalar(dimLength, Zero)
        )
    );
    volScalarField& cellDistance = cellDistancePtr_();

    // Internal field
    {
        const auto& cc = mesh.cellCentres();
        scalarField& fld = cellDistance.primitiveFieldRef();

        forAll(cc, i)
        {
            fld[i] = plane_.signedDistance(cc[i]);
        }
    }

    volScalarField::Boundary& cellDistanceBf =
        cellDistance.boundaryFieldRef();

    // Patch fields
    {
        forAll(cellDistanceBf, patchi)
        {
            if
            (
                isA<emptyFvPatchScalarField>
                (
                    cellDistanceBf[patchi]
                )
            )
            {
                cellDistanceBf.set
                (
                    patchi,
                    new calculatedFvPatchScalarField
                    (
                        mesh.boundary()[patchi],
                        cellDistance
                    )
                );

                const polyPatch& pp = mesh.boundary()[patchi].patch();
                pointField::subField cc = pp.patchSlice(mesh.faceCentres());

                fvPatchScalarField& fld = cellDistanceBf[patchi];
                fld.setSize(pp.size());
                forAll(fld, i)
                {
                    fld[i] = plane_.signedDistance(cc[i]);
                }
            }
            else
            {
                // Other side cell centres?
                const pointField& cc = mesh.C().boundaryField()[patchi];
                fvPatchScalarField& fld = cellDistanceBf[patchi];

                forAll(fld, i)
                {
                    fld[i] = plane_.signedDistance(cc[i]);
                }
            }
        }
    }


    // On processor patches the mesh.C() will already be the cell centre
    // on the opposite side so no need to swap cellDistance.


    // Distance to points
    pointDistance_.setSize(mesh.nPoints());
    {
        const pointField& pts = mesh.points();

        forAll(pointDistance_, i)
        {
            pointDistance_[i] = plane_.signedDistance(pts[i]);
        }
    }


    if (debug)
    {
        Pout<< "Writing cell distance:" << cellDistance.objectPath() << endl;
        cellDistance.write();
        pointScalarField pDist
        (
            IOobject
            (
                "pointDistance",
                mesh.time().timeName(),
                mesh.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            pointMesh::New(mesh),
            dimensionedScalar(dimLength, Zero)
        );
        pDist.primitiveFieldRef() = pointDistance_;

        Pout<< "Writing point distance:" << pDist.objectPath() << endl;
        pDist.write();
    }


    // Direct from cell field and point field.
    if (isoAlgo_ == isoSurfaceBase::ALGO_CELL)
    {
        isoSurfCellPtr_.reset
        (
            new isoSurfaceCell
            (
                fvm,
                cellDistance,
                pointDistance_,
                0,
                filter_,
                bounds_,
                mergeTol_
            )
        );
    }
    else if (isoAlgo_ == isoSurfaceBase::ALGO_TOPO)
    {
        isoSurfTopoPtr_.reset
        (
            new isoSurfaceTopo
            (
                fvm,
                cellDistance,
                pointDistance_,
                0,
                filter_,
                bounds_
            )
        );
    }
    else
    {
        isoSurfPtr_.reset
        (
            new isoSurface
            (
                cellDistance,
                pointDistance_,
                0,
                filter_,
                bounds_,
                mergeTol_
            )
        );
    }

    if (debug)
    {
        print(Pout);
        Pout<< endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sampledCuttingPlane::sampledCuttingPlane
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledSurface(name, mesh, dict),
    plane_(dict),
    bounds_(dict.getOrDefault("bounds", boundBox::invertedBox)),
    mergeTol_(dict.getOrDefault<scalar>("mergeTol", 1e-6)),
    isoAlgo_
    (
        isoSurfaceBase::algorithmNames.getOrDefault
        (
            "isoAlgorithm",
            dict,
            isoSurfaceBase::ALGO_POINT
        )
    ),
    filter_
    (
        isoSurfaceBase::getFilterType
        (
            dict,
            isoSurfaceBase::filterType::DIAGCELL
        )
    ),
    average_(dict.getOrDefault("average", false)),
    zoneNames_(),
    exposedPatchName_(),
    needsUpdate_(true),
    subMeshPtr_(nullptr),
    cellDistancePtr_(nullptr),
    isoSurfPtr_(nullptr),
    isoSurfCellPtr_(nullptr),
    isoSurfTopoPtr_(nullptr)
{
    if (!dict.readIfPresent("zones", zoneNames_) && dict.found("zone"))
    {
        zoneNames_.resize(1);
        dict.readEntry("zone", zoneNames_.first());
    }

    if (-1 != mesh.cellZones().findIndex(zoneNames_))
    {
        dict.readEntry("exposedPatchName", exposedPatchName_);

        if (-1 == mesh.boundaryMesh().findPatchID(exposedPatchName_))
        {
            FatalIOErrorInFunction(dict)
                << "Cannot find patch " << exposedPatchName_
                << " in which to put exposed faces." << endl
                << "Valid patches are " << mesh.boundaryMesh().names()
                << exit(FatalError);
        }

        DebugInfo
            << "Restricting to cellZone(s) " << flatOutput(zoneNames_)
            << " with exposed internal faces into patch "
            << exposedPatchName_ << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sampledCuttingPlane::needsUpdate() const
{
    return needsUpdate_;
}


bool Foam::sampledCuttingPlane::expire()
{
    if (debug)
    {
        Pout<< "sampledCuttingPlane::expire :"
            << " needsUpdate:" << needsUpdate_ << endl;
    }

    // Clear derived data
    clearGeom();

    // Already marked as expired
    if (needsUpdate_)
    {
        return false;
    }

    needsUpdate_ = true;
    return true;
}


bool Foam::sampledCuttingPlane::update()
{
    if (debug)
    {
        Pout<< "sampledCuttingPlane::update :"
            << " needsUpdate:" << needsUpdate_ << endl;
    }

    if (!needsUpdate_)
    {
        return false;
    }

    createGeometry();

    needsUpdate_ = false;
    return true;
}


Foam::tmp<Foam::scalarField>
Foam::sampledCuttingPlane::sample
(
    const interpolation<scalar>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::vectorField>
Foam::sampledCuttingPlane::sample
(
    const interpolation<vector>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::sphericalTensorField>
Foam::sampledCuttingPlane::sample
(
    const interpolation<sphericalTensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::symmTensorField>
Foam::sampledCuttingPlane::sample
(
    const interpolation<symmTensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::tensorField>
Foam::sampledCuttingPlane::sample
(
    const interpolation<tensor>& sampler
) const
{
    return sampleOnFaces(sampler);
}


Foam::tmp<Foam::scalarField>
Foam::sampledCuttingPlane::interpolate
(
    const interpolation<scalar>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::vectorField>
Foam::sampledCuttingPlane::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::sphericalTensorField>
Foam::sampledCuttingPlane::interpolate
(
    const interpolation<sphericalTensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::symmTensorField>
Foam::sampledCuttingPlane::interpolate
(
    const interpolation<symmTensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


Foam::tmp<Foam::tensorField>
Foam::sampledCuttingPlane::interpolate
(
    const interpolation<tensor>& interpolator
) const
{
    return sampleOnPoints(interpolator);
}


void Foam::sampledCuttingPlane::print(Ostream& os) const
{
    os  << "sampledCuttingPlane: " << name() << " :"
        << "  plane:" << plane_
        << "  faces:" << faces().size()
        << "  points:" << points().size();
}


// ************************************************************************* //
