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

#include "distanceSurface.H"
#include "dictionary.H"
#include "volFields.H"
#include "volPointInterpolation.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(distanceSurface, 0);
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
    static isoSurfaceBase::algorithmType getIsoAlgorithm(const dictionary& dict)
    {
        // Previously 'cell' (bool), but now 'isoAlgorithm' (enum)

        // Default (bool) for 1906 and earlier
        bool useCell = true;

        // Default (enum) after 1906
        isoSurfaceBase::algorithmType algo = isoSurfaceBase::ALGO_CELL;

        if
        (
            !isoSurfaceBase::algorithmNames.readIfPresent
            (
                "isoAlgorithm", dict, algo
            )
            // When above fails, use 'compat' to also get upgrade messages
         && dict.readIfPresentCompat
            (
                "isoAlgorithm", {{"cell", 1906}}, useCell
            )
        )
        {
            return
            (
                useCell
              ? isoSurfaceBase::ALGO_CELL
              : isoSurfaceBase::ALGO_POINT
            );
        }

        return algo;
    }

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distanceSurface::distanceSurface
(
    const word& defaultSurfaceName,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    surfPtr_
    (
        searchableSurface::New
        (
            dict.get<word>("surfaceType"),
            IOobject
            (
                dict.getOrDefault("surfaceName", defaultSurfaceName),
                mesh.time().constant(), // directory
                "triSurface",           // instance
                mesh.time(),            // registry
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            dict
        )
    ),
    distance_(dict.get<scalar>("distance")),
    signed_
    (
        distance_ < 0 || equal(distance_, Zero) || dict.get<bool>("signed")
    ),
    isoAlgo_(getIsoAlgorithm(dict)),
    filter_
    (
        isoSurfaceBase::getFilterType
        (
            dict,
            isoSurfaceBase::filterType::DIAGCELL
        )
    ),
    bounds_(dict.getOrDefault("bounds", boundBox::invertedBox)),
    isoSurfPtr_(nullptr),
    isoSurfCellPtr_(nullptr),
    isoSurfTopoPtr_(nullptr)
{}


Foam::distanceSurface::distanceSurface
(
    const polyMesh& mesh,
    const bool interpolate,
    const word& surfaceType,
    const word& surfaceName,
    const scalar distance,
    const bool signedDistance,
    const isoSurfaceBase::algorithmType algo,
    const isoSurfaceBase::filterType filter,
    const boundBox& bounds
)
:
    mesh_(mesh),
    surfPtr_
    (
        searchableSurface::New
        (
            surfaceType,
            IOobject
            (
                surfaceName,            // name
                mesh.time().constant(), // directory
                "triSurface",           // instance
                mesh.time(),            // registry
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            dictionary()
        )
    ),
    distance_(distance),
    signed_
    (
        signedDistance || distance_ < 0 || equal(distance_, Zero)
    ),
    isoAlgo_(algo),
    filter_(filter),
    bounds_(bounds),
    isoSurfPtr_(nullptr),
    isoSurfCellPtr_(nullptr),
    isoSurfTopoPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::distanceSurface::createGeometry()
{
    if (debug)
    {
        Pout<< "distanceSurface::createGeometry updating geometry." << endl;
    }

    // Clear any stored topologies
    isoSurfPtr_.clear();
    isoSurfCellPtr_.clear();
    isoSurfTopoPtr_.clear();

    const fvMesh& fvm = static_cast<const fvMesh&>(mesh_);

    // Distance to cell centres
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    cellDistancePtr_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "distanceSurface.cellDistance",
                fvm.time().timeName(),
                fvm.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            fvm,
            dimensionedScalar(dimLength, Zero)
        )
    );
    volScalarField& cellDistance = *cellDistancePtr_;

    // For distance = 0 (and isoSurfaceCell) we apply additional filtering
    // to limit the extent of open edges.

    const bool isZeroDist  = equal(distance_, Zero);
    const bool filterCells =
    (
        isZeroDist
     && isoAlgo_ != isoSurfaceBase::ALGO_POINT
    );

    bitSet ignoreCells;
    if (filterCells)
    {
        ignoreCells.resize(fvm.nCells());
    }

    // Internal field
    {
        const pointField& cc = fvm.C();
        scalarField& fld = cellDistance.primitiveFieldRef();

        List<pointIndexHit> nearest;
        surfPtr_().findNearest
        (
            cc,
            scalarField(cc.size(), GREAT),
            nearest
        );

        if (signed_ || isZeroDist)
        {
            vectorField norms;
            surfPtr_().getNormal(nearest, norms);

            boundBox cellBb;

            forAll(norms, i)
            {
                const point diff(cc[i] - nearest[i].hitPoint());

                fld[i] =
                (
                    isZeroDist // Use normal distance
                  ? (diff & norms[i])
                  : Foam::sign(diff & norms[i]) * Foam::mag(diff)
                );

                if (filterCells)
                {
                    cellBb.clear();
                    cellBb.add(fvm.points(), fvm.cellPoints(i));

                    // Expand slightly to catch corners
                    cellBb.inflate(0.1);

                    if (!cellBb.contains(nearest[i].hitPoint()))
                    {
                        ignoreCells.set(i);
                    }
                }
            }
        }
        else
        {
            forAll(nearest, i)
            {
                fld[i] = Foam::mag(cc[i] - nearest[i].hitPoint());

                // No filtering for unsigned or distance != 0.
            }
        }
    }

    // Patch fields
    {
        volScalarField::Boundary& cellDistanceBf =
            cellDistance.boundaryFieldRef();

        forAll(fvm.C().boundaryField(), patchi)
        {
            const pointField& cc = fvm.C().boundaryField()[patchi];
            fvPatchScalarField& fld = cellDistanceBf[patchi];

            List<pointIndexHit> nearest;
            surfPtr_().findNearest
            (
                cc,
                scalarField(cc.size(), GREAT),
                nearest
            );

            if (signed_)
            {
                vectorField norms;
                surfPtr_().getNormal(nearest, norms);

                forAll(norms, i)
                {
                    const point diff(cc[i] - nearest[i].hitPoint());

                    fld[i] =
                    (
                        isZeroDist // Use normal distance
                      ? (diff & norms[i])
                      : Foam::sign(diff & norms[i]) * Foam::mag(diff)
                    );
                }
            }
            else
            {
                forAll(nearest, i)
                {
                    fld[i] = Foam::mag(cc[i] - nearest[i].hitPoint());
                }
            }
        }
    }


    // On processor patches the mesh.C() will already be the cell centre
    // on the opposite side so no need to swap cellDistance.


    // Distance to points
    pointDistance_.setSize(fvm.nPoints());
    {
        const pointField& pts = fvm.points();

        List<pointIndexHit> nearest;
        surfPtr_().findNearest
        (
            pts,
            scalarField(pts.size(), GREAT),
            nearest
        );

        if (signed_)
        {
            vectorField norms;
            surfPtr_().getNormal(nearest, norms);

            forAll(norms, i)
            {
                const point diff(pts[i] - nearest[i].hitPoint());

                pointDistance_[i] =
                (
                    isZeroDist // Use normal distance
                  ? (diff & norms[i])
                  : Foam::sign(diff & norms[i]) * Foam::mag(diff)
                );
            }
        }
        else
        {
            forAll(nearest, i)
            {
                pointDistance_[i] = Foam::mag(pts[i] - nearest[i].hitPoint());
            }
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
                "distanceSurface.pointDistance",
                fvm.time().timeName(),
                fvm.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            pointMesh::New(fvm),
            dimensionedScalar(dimLength, Zero)
        );
        pDist.primitiveFieldRef() = pointDistance_;

        Pout<< "Writing point distance:" << pDist.objectPath() << endl;
        pDist.write();
    }

    // Don't need ignoreCells if there is nothing to ignore.
    if (!ignoreCells.any())
    {
        ignoreCells.clear();
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
                distance_,
                filter_,
                bounds_,
                1e-6,  // mergeTol
                ignoreCells
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
                distance_,
                filter_,
                bounds_,
                ignoreCells
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
                distance_,
                filter_,
                bounds_,
                1e-6
            )
        );
    }

    if (debug)
    {
        print(Pout);
        Pout<< endl;
    }
}


void Foam::distanceSurface::print(Ostream& os) const
{
    os  << "  surface:" << surfaceName()
        << "  distance:" << distance()
        << "  faces:" << surface().surfFaces().size()
        << "  points:" << surface().points().size();
}


// ************************************************************************* //
