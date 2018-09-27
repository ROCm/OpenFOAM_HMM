/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017-2018 OpenCFD Ltd.
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

#include "faceShading.H"
#include "fvMesh.H"
#include "boundaryRadiationProperties.H"
#include "OFstream.H"
#include "cyclicAMIPolyPatch.H"
#include "volFields.H"
#include "distributedTriSurfaceMesh.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faceShading, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::faceShading::writeRays
(
    const fileName& fName,
    const DynamicField<point>& endCf,
    const pointField& myFc
)
{
    OFstream str(fName);
    label vertI = 0;

    Pout<< "Dumping rays to " << str.name() << endl;

    forAll(myFc, faceI)
    {
        meshTools::writeOBJ(str, myFc[faceI]);
        vertI++;
        meshTools::writeOBJ(str, endCf[faceI]);
        vertI++;
        str << "l " << vertI-1 << ' ' << vertI << nl;
    }
    str.flush();

    Pout<< "cmd: objToVTK " << fName.c_str() << endl;

    stringList cmd({"objToVTK", fName, fName.lessExt().ext("vtk")});
    Foam::system(cmd);
}


Foam::triSurface Foam::faceShading::triangulate
(
    const labelHashSet& includePatches,
    const List<labelHashSet>& includeAllFacesPerPatch
)
{
    const polyBoundaryMesh& bMesh = mesh_.boundaryMesh();

    // Storage for surfaceMesh. Size estimate.
    DynamicList<labelledTri> triangles(mesh_.nBoundaryFaces());

    label newPatchI = 0;

    for (const label patchI : includePatches)
    {
        const polyPatch& patch = bMesh[patchI];
        const pointField& points = patch.points();

        label nTriTotal = 0;

        if (includeAllFacesPerPatch[patchI].size())
        {
            for (const label patchFaceI : includeAllFacesPerPatch[patchI])
            {
                const face& f = patch[patchFaceI];

                faceList triFaces(f.nTriangles(points));

                label nTri = 0;

                f.triangles(points, nTri, triFaces);

                forAll(triFaces, triFaceI)
                {
                    const face& f = triFaces[triFaceI];

                    triangles.append
                    (
                        labelledTri(f[0], f[1], f[2], newPatchI)
                    );
                    nTriTotal++;
                }
            }
            newPatchI++;
        }
    }

    triangles.shrink();

    // Create globally numbered tri surface
    triSurface rawSurface(triangles, mesh_.points());

    // Create locally numbered tri surface
    triSurface surface
    (
        rawSurface.localFaces(),
        rawSurface.localPoints()
    );

    // Add patch names to surface
    surface.patches().setSize(newPatchI);

    newPatchI = 0;

    for (const label patchI : includePatches)
    {
        const polyPatch& patch = bMesh[patchI];

        if (includeAllFacesPerPatch[patchI].size())
        {
            surface.patches()[newPatchI].name() = patch.name();
            surface.patches()[newPatchI].geometricType() = patch.type();

            newPatchI++;
        }
    }

    return surface;
}


void Foam::faceShading::calculate()
{
    const radiation::boundaryRadiationProperties& boundaryRadiation =
        radiation::boundaryRadiationProperties::New(mesh_);

    label nFaces = 0;          //total number of direct hit faces

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    DynamicList<point> dynCf(nFaces);
    DynamicList<label> dynFacesI;

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        const pointField& cf = pp.faceCentres();

        if (!pp.coupled() && !isA<cyclicAMIPolyPatch>(pp))
        {
            const tmp<scalarField> tt =
                boundaryRadiation.transmissivity(patchI);
            const scalarField& t = tt();
            const vectorField& n = pp.faceNormals();

            forAll(n, faceI)
            {
                const vector nf(n[faceI]);
                if (((direction_ & nf) > 0) && (t[faceI] == 0.0))
                {
                    dynFacesI.append(faceI + pp.start());
                    dynCf.append(cf[faceI]);
                    nFaces++;
                }
            }
        }
    }

    label numberPotentialHits = nFaces;

    reduce(numberPotentialHits, sumOp<label>());

    Info<< "Number of 'potential' direct hits : "
        << numberPotentialHits << endl;

    labelList hitFacesIds(nFaces);
    hitFacesIds.transfer(dynFacesI);

    pointField Cfs(hitFacesIds.size());
    Cfs.transfer(dynCf);

    // * * * * * * * * * * * * * * *
    // Create distributedTriSurfaceMesh
    Random rndGen(653213);

    // Determine mesh bounding boxes:
    List<treeBoundBox> meshBb
    (
        1,
        treeBoundBox
        (
            boundBox(mesh_.points(), false)
        ).extend(rndGen, 1e-3)
    );

    // Dummy bounds dictionary
    dictionary dict;
    dict.add("bounds", meshBb);
    dict.add
    (
        "distributionType",
        distributedTriSurfaceMesh::distributionTypeNames_
        [
            distributedTriSurfaceMesh::FROZEN
        ]
    );
    dict.add("mergeDistance", SMALL);

    labelHashSet includePatches;
    List<labelHashSet> includeAllFacesPerPatch(patches.size());

    forAll(patches, patchI)
    {
        const polyPatch& pp = patches[patchI];
        if (!pp.coupled() && !isA<cyclicAMIPolyPatch>(pp))
        {
            includePatches.insert(patchI);

            const tmp<scalarField> tt =
                boundaryRadiation.transmissivity(patchI);
            const scalarField& tau = tt();

            forAll(pp, faceI)
            {
                if (tau[faceI] == 0.0)
                {
                    includeAllFacesPerPatch[patchI].insert
                    (
                        faceI
                    );
                }
            }
        }
    }

    triSurface localSurface = triangulate
    (
        includePatches,
        includeAllFacesPerPatch
    );

    distributedTriSurfaceMesh surfacesMesh
    (
        IOobject
        (
            "opaqueSurface.stl",
            mesh_.time().constant(),    // directory
            "triSurface",               // instance
            mesh_.time(),               // registry
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        localSurface,
        dict
    );

    if (debug)
    {
        surfacesMesh.searchableSurface::write();
    }

    scalar maxBounding = 5.0*mag(mesh_.bounds().max() - mesh_.bounds().min());

    reduce(maxBounding, maxOp<scalar>());

    // Calculate index of faces which have a direct hit (local)
    DynamicList<label> rayStartFace(nFaces + 0.01*nFaces);

    // Shoot Rays
    // * * * * * * * * * * * * * * * *
    {

        DynamicField<point> start(nFaces);
        DynamicField<point> end(start.size());
        DynamicList<label> startIndex(start.size());

        label i = 0;
        do
        {
            for (; i < Cfs.size(); i++)
            {
                const point& fc = Cfs[i];

                const label myFaceId = hitFacesIds[i];

                const vector d(direction_*maxBounding);

                start.append(fc - 0.001*d);

                startIndex.append(myFaceId);

                end.append(fc - d);

            }

        }while (returnReduce(i < Cfs.size(), orOp<bool>()));

        List<pointIndexHit> hitInfo(startIndex.size());
        surfacesMesh.findLine(start, end, hitInfo);

        // Collect the rays which has 'only one not wall' obstacle between
        // start and end.
        // If the ray hit itself get stored in dRayIs
        forAll(hitInfo, rayI)
        {
            if (!hitInfo[rayI].hit())
            {
                rayStartFace.append(startIndex[rayI]);
            }
        }

        // Plot all rays between visible faces.
        if (debug)
        {
            writeRays
            (
                mesh_.time().path()/"allVisibleFaces.obj",
                end,
                start
            );
        }

        start.clear();
        startIndex.clear();
        end.clear();
    }

    rayStartFaces_.transfer(rayStartFace);

    if (debug)
    {
        tmp<volScalarField> thitFaces
        (
            new volScalarField
            (
                IOobject
                (
                    "hitFaces",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar(dimless, Zero)
            )
        );

        volScalarField& hitFaces = thitFaces.ref();
        volScalarField::Boundary& hitFacesBf = hitFaces.boundaryFieldRef();

        hitFacesBf = 0.0;
        forAll(rayStartFaces_, i)
        {
            const label faceI = rayStartFaces_[i];
            label patchID = patches.whichPatch(faceI);
            const polyPatch& pp = patches[patchID];
            hitFacesBf[patchID][faceI - pp.start()] = 1.0;
        }
        hitFaces.write();
    }

    label totalHitFaces = rayStartFaces_.size();

    reduce(totalHitFaces, sumOp<label>());

    Info<< "Total number of hit faces : " <<  totalHitFaces << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faceShading::faceShading
(
    const fvMesh& mesh,
    const vector dir,
    const labelList& hitFaceList
)
:
    mesh_(mesh),
    direction_(dir),
    rayStartFaces_(hitFaceList)
{}



Foam::faceShading::faceShading
(
    const fvMesh& mesh,
    const vector dir
)
:
    mesh_(mesh),
    direction_(dir),
    rayStartFaces_(0)
{
    calculate();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::faceShading::~faceShading()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::faceShading::correct()
{
    calculate();
}


// ************************************************************************* //
