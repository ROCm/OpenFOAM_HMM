/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "distanceSurface.H"
#include "dictionary.H"
#include "volFields.H"
#include "volPointInterpolation.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
#include "isoSurfaceCell.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(distanceSurface, 0);
    addToRunTimeSelectionTable(sampledSurface, distanceSurface, word);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::distanceSurface::createGeometry() const
{
    // Clear any stored topo
    facesPtr_.clear();

    // Distance to cell centres
    scalarField cellDistance(mesh().nCells());
    {
        List<pointIndexHit> nearest;
        surfPtr_().findNearest
        (
            mesh().cellCentres(),
            scalarField(mesh().nCells(), GREAT),
            nearest
        );

        if (signed_)
        {
            vectorField normal;
            surfPtr_().getNormal(nearest, normal);

            forAll(cellDistance, cellI)
            {
                vector d(mesh().cellCentres()[cellI]-nearest[cellI].hitPoint());

                if ((d&normal[cellI]) > 0)
                {
                    cellDistance[cellI] = Foam::mag(d);
                }
                else
                {
                    cellDistance[cellI] = -Foam::mag(d);
                }
            }
        }
        else
        {
            forAll(cellDistance, cellI)
            {
                cellDistance[cellI] = Foam::mag
                (
                    nearest[cellI].hitPoint()
                  - mesh().cellCentres()[cellI]
                );
            }
        }
    }

    // Distance to points
    scalarField pointDistance(mesh().nPoints());
    {
        List<pointIndexHit> nearest;
        surfPtr_().findNearest
        (
            mesh().points(),
            scalarField(mesh().nPoints(), GREAT),
            nearest
        );
        forAll(pointDistance, pointI)
        {
            pointDistance[pointI] = Foam::mag
            (
                nearest[pointI].hitPoint()
              - mesh().points()[pointI]
            );
        }
    }

    //- Direct from cell field and point field.
    const isoSurfaceCell iso
    (
        mesh(),
        cellDistance,
        pointDistance,
        distance_,
        regularise_
    );

    ////- From point field and interpolated cell.
    //scalarField cellAvg(mesh().nCells(), scalar(0.0));
    //labelField nPointCells(mesh().nCells(), 0);
    //{
    //    for (label pointI = 0; pointI < mesh().nPoints(); pointI++)
    //    {
    //        const labelList& pCells = mesh().pointCells(pointI);
    //
    //        forAll(pCells, i)
    //        {
    //            label cellI = pCells[i];
    //
    //            cellAvg[cellI] += pointDistance[pointI];
    //            nPointCells[cellI]++;
    //        }
    //    }
    //}
    //forAll(cellAvg, cellI)
    //{
    //    cellAvg[cellI] /= nPointCells[cellI];
    //}
    //
    //const isoSurface iso
    //(
    //    mesh(),
    //    cellAvg,
    //    pointDistance,
    //    distance_,
    //    regularise_
    //);


    const_cast<distanceSurface&>(*this).triSurface::operator=(iso);
    meshCells_ = iso.meshCells();

    if (debug)
    {
        print(Pout);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distanceSurface::distanceSurface
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    sampledSurface(name, mesh, dict),
    surfPtr_
    (
        searchableSurface::New
        (
            dict.lookup("surfaceType"),
            IOobject
            (
                dict.lookupOrDefault("surfaceName", name),  // name
                mesh.time().constant(),                     // directory
                "triSurface",                               // instance
                mesh.time(),                                // registry
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            dict
        )
    ),
    distance_(readScalar(dict.lookup("distance"))),
    signed_(readBool(dict.lookup("signed"))),
    regularise_(dict.lookupOrDefault("regularise", true)),
    zoneName_(word::null),
    facesPtr_(NULL),
    storedTimeIndex_(-1),
    meshCells_(0)
{
//    label zoneId = -1;
//    if (dict.readIfPresent("zone", zoneName_))
//    {
//        zoneId = mesh.cellZones().findZoneID(zoneName_);
//        if (debug && zoneId < 0)
//        {
//            Info<< "cellZone \"" << zoneName_
//                << "\" not found - using entire mesh"
//                << endl;
//        }
//    }
    correct(true);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::distanceSurface::~distanceSurface()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::distanceSurface::correct(const bool meshChanged)
{
    // Only change of mesh changes plane - zone restriction gets lost
    if (meshChanged)
    {
        createGeometry();
    }
}


Foam::tmp<Foam::scalarField>
Foam::distanceSurface::sample
(
    const volScalarField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::vectorField>
Foam::distanceSurface::sample
(
    const volVectorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::sphericalTensorField>
Foam::distanceSurface::sample
(
    const volSphericalTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::symmTensorField>
Foam::distanceSurface::sample
(
    const volSymmTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::tensorField>
Foam::distanceSurface::sample
(
    const volTensorField& vField
) const
{
    return sampleField(vField);
}


Foam::tmp<Foam::scalarField>
Foam::distanceSurface::interpolate
(
    const interpolation<scalar>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::vectorField>
Foam::distanceSurface::interpolate
(
    const interpolation<vector>& interpolator
) const
{
    return interpolateField(interpolator);
}

Foam::tmp<Foam::sphericalTensorField>
Foam::distanceSurface::interpolate
(
    const interpolation<sphericalTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::symmTensorField>
Foam::distanceSurface::interpolate
(
    const interpolation<symmTensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


Foam::tmp<Foam::tensorField>
Foam::distanceSurface::interpolate
(
    const interpolation<tensor>& interpolator
) const
{
    return interpolateField(interpolator);
}


void Foam::distanceSurface::print(Ostream& os) const
{
    os  << "distanceSurface: " << name() << " :"
        << "  surface:" << surfPtr_().name()
        << "  distance:" << distance_
        << "  faces:" << faces().size()
        << "  points:" << points().size();
}


// ************************************************************************* //
