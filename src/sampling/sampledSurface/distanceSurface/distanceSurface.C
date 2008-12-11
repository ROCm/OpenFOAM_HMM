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
#include "isoSurface.H"
//#include "isoSurfaceCell.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(distanceSurface, 0);
    addToRunTimeSelectionTable(sampledSurface, distanceSurface, word);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::distanceSurface::createGeometry()
{
    // Clear any stored topo
    facesPtr_.clear();


    const fvMesh& fvm = static_cast<const fvMesh&>(mesh());

    // Distance to cell centres
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    cellDistancePtr_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "cellDistance",
                fvm.time().timeName(),
                fvm.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            fvm,
            dimensionedScalar("zero", dimless/dimTime, 0)
        )
    );
    volScalarField& cellDistance = cellDistancePtr_();

    // Internal field
    {
        const pointField& cc = fvm.C();
        scalarField& fld = cellDistance.internalField();

        List<pointIndexHit> nearest;
        surfPtr_().findNearest
        (
            cc,
            scalarField(cc.size(), GREAT),
            nearest
        );

        if (signed_)
        {
            vectorField normal;
            surfPtr_().getNormal(nearest, normal);

            forAll(nearest, i)
            {
                vector d(cc[i]-nearest[i].hitPoint());

                if ((d&normal[i]) > 0)
                {
                    fld[i] = Foam::mag(d);
                }
                else
                {
                    fld[i] = -Foam::mag(d);
                }
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

    // Patch fields
    {
        forAll(fvm.C().boundaryField(), patchI)
        {
            const pointField& cc = fvm.C().boundaryField()[patchI];
            fvPatchScalarField& fld = cellDistance.boundaryField()[patchI];

            List<pointIndexHit> nearest;
            surfPtr_().findNearest
            (
                cc,
                scalarField(cc.size(), GREAT),
                nearest
            );

            if (signed_)
            {
                vectorField normal;
                surfPtr_().getNormal(nearest, normal);

                forAll(nearest, i)
                {
                    vector d(cc[i]-nearest[i].hitPoint());

                    if ((d&normal[i]) > 0)
                    {
                        fld[i] = Foam::mag(d);
                    }
                    else
                    {
                        fld[i] = -Foam::mag(d);
                    }
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
        List<pointIndexHit> nearest;
        surfPtr_().findNearest
        (
            fvm.points(),
            scalarField(fvm.nPoints(), GREAT),
            nearest
        );
        forAll(pointDistance_, pointI)
        {
            pointDistance_[pointI] = Foam::mag
            (
                nearest[pointI].hitPoint()
              - fvm.points()[pointI]
            );
        }
    }

    //- Direct from cell field and point field.
    isoSurfPtr_.reset
    (
        new isoSurface
        (
            cellDistance,
            pointDistance_,
            distance_,
            regularise_
        )
    );

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
    isoSurfPtr_(NULL),
    facesPtr_(NULL)
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
