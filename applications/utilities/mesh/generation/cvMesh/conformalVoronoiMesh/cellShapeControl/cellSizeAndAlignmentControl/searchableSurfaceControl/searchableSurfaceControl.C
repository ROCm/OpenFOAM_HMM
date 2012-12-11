/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

#include "searchableSurfaceControl.H"
#include "addToRunTimeSelectionTable.H"
#include "cellSizeFunction.H"
#include "triSurfaceMesh.H"
#include "searchableBox.H"
#include "tetrahedron.H"
#include "vectorTools.H"
#include "quaternion.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(searchableSurfaceControl, 0);
addToRunTimeSelectionTable
(
    cellSizeAndAlignmentControl,
    searchableSurfaceControl,
    dictionary
);

}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

//Foam::tensor Foam::surfaceControl::requiredAlignment
//(
//    const Foam::point& pt,
//    const vectorField& ptNormals
//) const
//{
////    pointIndexHit surfHit;
////    label hitSurface;
////
////    geometryToConformTo_.findSurfaceNearest
////    (
////        pt,
////        sqr(GREAT),
////        surfHit,
////        hitSurface
////    );
////
////    if (!surfHit.hit())
////    {
////        FatalErrorIn
////        (
////            "Foam::tensor Foam::conformalVoronoiMesh::requiredAlignment"
////        )
////            << "findSurfaceNearest did not find a hit across the surfaces."
////            << exit(FatalError) << endl;
////    }
////
//////     Primary alignment
////
////    vectorField norm(1);
////
////    allGeometry_[hitSurface].getNormal
////    (
////        List<pointIndexHit>(1, surfHit),
////        norm
////    );
////
////    const vector np = norm[0];
////
////    const tensor Rp = rotationTensor(vector(0,0,1), np);
////
////    return (Rp);
//
////    Info<< "Point : " << pt << endl;
////    forAll(ptNormals, pnI)
////    {
////        Info<< "    normal " << pnI << " : " << ptNormals[pnI] << endl;
////    }
//
//    vector np = ptNormals[0];
//
//    const tensor Rp = rotationTensor(vector(0,0,1), np);
//
//    vector na = vector::zero;
//
//    scalar smallestAngle = GREAT;
//
//    for (label pnI = 1; pnI < ptNormals.size(); ++pnI)
//    {
//        const vector& nextNormal = ptNormals[pnI];
//
//        const scalar cosPhi = vectorTools::cosPhi(np, nextNormal);
//
//        if (mag(cosPhi) < smallestAngle)
//        {
//            na = nextNormal;
//            smallestAngle = cosPhi;
//        }
//    }
//
//    // Secondary alignment
//    vector ns = np ^ na;
//
//    if (mag(ns) < SMALL)
//    {
//        WarningIn("conformalVoronoiMesh::requiredAlignment")
//            << "Parallel normals detected in spoke search." << nl
//            << "point: " << pt << nl
//            << "np   : " << np << nl
//            << "na   : " << na << nl
//            << "ns   : " << ns << nl
//            << endl;
//
//        ns = np;
//    }
//
//    ns /= mag(ns);
//
//    tensor Rs = rotationTensor((Rp & vector(0,1,0)), ns);
//
////    Info<< "Point " << pt << nl
////        << "      np : " << np << nl
////        << "      ns : " << ns << nl
////        << "      Rp : " << Rp << nl
////        << "      Rs : " << Rs << nl
////        << "    Rs&Rp: " << (Rs & Rp) << endl;
//
//    return (Rs & Rp);
//}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::searchableSurfaceControl::searchableSurfaceControl
(
    const Time& runTime,
    const word& name,
    const dictionary& controlFunctionDict,
    const conformationSurfaces& allGeometry
)
:
    cellSizeAndAlignmentControl(runTime, name, controlFunctionDict, allGeometry),
    searchableSurface_(allGeometry.geometry()[name]),
    allGeometry_(allGeometry),
    cellSizeFunction_
    (
        cellSizeFunction::New(controlFunctionDict, searchableSurface_)
    )
//    geometryToConformTo_(geometryToConformTo),
//    surfaces_(),
//    cellSizeFunctions_(),
//    triangulatedMesh_()
{
//    const dictionary& surfacesDict = coeffDict();
//
//    Info<< nl << "Reading cellSizeControlGeometry" << endl;
//
//    surfaces_.setSize(surfacesDict.size());
//
//    cellSizeFunctions_.setSize(surfacesDict.size());
//
//    label surfI = 0;
//
//    DynamicList<point> pointsToInsert;
//    DynamicList<scalar> sizesToInsert;
//    DynamicList<tensor> alignmentsToInsert;
//
//    forAllConstIter(dictionary, surfacesDict, iter)
//    {
//        const dictionary& surfaceSubDict
//        (
//            surfacesDict.subDict(iter().keyword())
//        );
//
//        // If the "surface" keyword is not found in the dictionary, assume that
//        // the name of the dictionary is the surface. Distinction required to
//        // allow the same surface to be used multiple times to supply multiple
//        // cellSizeFunctions
//
//        word surfaceName = surfaceSubDict.lookupOrDefault<word>
//        (
//            "surface",
//            iter().keyword()
//        );
//
//        surfaces_[surfI] = allGeometry_.findSurfaceID(surfaceName);
//
//        if (surfaces_[surfI] < 0)
//        {
//            FatalErrorIn
//            (
//                "Foam::surfaceControl::surfaceControl"
//            )   << "No surface " << surfaceName << " found. "
//                << "Valid geometry is " << nl << allGeometry_.names()
//                << exit(FatalError);
//        }
//
//        const searchableSurface& surface = allGeometry_[surfaces_[surfI]];
//
//        Info<< nl << "    " << iter().keyword() << nl
//            << "    surface: " << surfaceName << nl
//            << "    size   : " << surface.size() << endl;
//
//        cellSizeFunctions_.set
//        (
//            surfI,
//            cellSizeFunction::New
//            (
//                surfaceSubDict,
//                surface
//            )
//        );
//
//        surfI++;
//
//        if (isA<triSurfaceMesh>(surface))
//        {
//            const triSurfaceMesh& tsm
//                = refCast<const triSurfaceMesh>(surface);
//
//            const pointField& points = tsm.points();
//            const vectorField& faceNormals = tsm.faceNormals();
//            const labelListList& pointFaces = tsm.pointFaces();
//
//            Info<< "    Number of points: " << tsm.nPoints() << endl;
//
//            forAll(points, pI)
//            {
//                const Foam::point& pt = points[pI];
//                const labelList& ptFaces = pointFaces[pI];
//
//                vectorField pointNormals(ptFaces.size());
//
//                forAll(pointNormals, pnI)
//                {
//                    pointNormals[pnI] = faceNormals[ptFaces[pnI]];
//                }
//
//                pointsToInsert.append(pt);
//
//                // Get the value of the point from surfaceCellSizeFunction. If
//                // adding points internally then will need to interpolate.
//                scalar newSize = 0;
//
//                cellSizeFunctions_[surfI - 1].cellSize(pt, newSize);
//                sizesToInsert.append(newSize);
//
//                tensor newAlignment = requiredAlignment(pt, pointNormals);
//
//                alignmentsToInsert.append(newAlignment);
//            }
//        }
//    }
//
//    // Add the global bound box to ensure that all internal point queries
//    // will return sizes and alignments
////    boundBox bb = allGeometry_.bounds();
////
////    pointField bbPoints = bb.points();
////
////    forAll(bbPoints, pI)
////    {
////        pointsToInsert.append(bbPoints[pI]);
////        sizesToInsert.append(defaultCellSize());
////        alignmentsToInsert.append(tensor(1,0,0,0,1,0,0,0,1));
////    }
//
//    triangulatedMesh_.set
//    (
//        new triangulatedMesh
//        (
//            runTime,
//            pointsToInsert,
//            sizesToInsert,
//            alignmentsToInsert,
//            defaultCellSize()
//        )
//    );
//
//    scalar factor = 1.0;
//    label maxIteration = cellShapeControlDict.lookupOrDefault<label>
//    (
//        "maxRefinementIterations", 1
//    );
//
//
//    for (label iteration = 0; iteration < maxIteration; ++iteration)
//    {
//        Info<< "Iteration : " << iteration << endl;
//
//        label nRefined = triangulatedMesh_().refineTriangulation
//        (
//            factor,
//            allGeometry_,
//            geometryToConformTo_
//        );
//
//        Info<< "    Number of cells refined in refinement iteration : "
//            << nRefined << nl << endl;
//
//        if (nRefined <= 0 && iteration != 0)
//        {
//            break;
//        }
//
//        factor *= 1.5;
//    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::searchableSurfaceControl::~searchableSurfaceControl()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

//Foam::scalar Foam::surfaceControl::cellSize(const point& pt) const
//{
//    scalarList bary;
//    Cell_handle ch;
//
//    triangulatedMesh_().barycentricCoords(pt, bary, ch);
//
//    scalar size = 0;
//    forAll(bary, pI)
//    {
//        size += bary[pI]*ch->vertex(pI)->size();
//    }
//
//    return size;
//}
//
//
////- Return the cell alignment at the given location
//Foam::tensor Foam::surfaceControl::cellAlignment(const point& pt) const
//{
//    scalarList bary;
//    Cell_handle ch;
//
//    triangulatedMesh_().barycentricCoords(pt, bary, ch);
//
////    vectorField cartesianDirections(3);
////
////    cartesianDirections[0] = vector(0,0,1);
////    cartesianDirections[1] = vector(0,1,0);
////    cartesianDirections[2] = vector(1,0,0);
////
////    // Rearrange each alignment tensor so that the x/y/z components are
////    // in order of whichever makes the smallest angle with the global coordinate
////    // system
////    FixedList<tensor, 4> alignments;
////
////    forAll(alignments, aI)
////    {
////        tensor a = ch->vertex(aI)->alignment();
////
////        tensor tmpA = a;
////
//////        Info<< nl << indent<< a << endl;
////
////        scalar minAngle = 0;
////
////        scalar axx = vectorTools::cosPhi(a.x(), cartesianDirections[0]);
////        scalar axy = vectorTools::cosPhi(a.y(), cartesianDirections[0]);
////        scalar axz = vectorTools::cosPhi(a.z(), cartesianDirections[0]);
////
////        scalar ayx = vectorTools::cosPhi(a.x(), cartesianDirections[1]);
////        scalar ayy = vectorTools::cosPhi(a.y(), cartesianDirections[1]);
////        scalar ayz = vectorTools::cosPhi(a.z(), cartesianDirections[1]);
////
////        scalar azx = vectorTools::cosPhi(a.x(), cartesianDirections[2]);
////        scalar azy = vectorTools::cosPhi(a.y(), cartesianDirections[2]);
////        scalar azz = vectorTools::cosPhi(a.z(), cartesianDirections[2]);
////
//////        Info<< indent << axx << " " << axy << " " << axz << nl
//////            << indent << ayx << " " << ayy << " " << ayz << nl
//////            << indent << azx << " " << azy << " " << azz << endl;
////
////        if (mag(axx) >= minAngle)
////        {
////            tmpA.xx() = mag(a.xx()); tmpA.xy() = mag(a.xy()); tmpA.xz() = mag(a.xz());
////            minAngle = mag(axx);
////        }
////        if (mag(axy) >= minAngle)
////        {
////            tmpA.xx() = mag(a.yx()); tmpA.xy() = mag(a.yy()); tmpA.xz() = mag(a.yz());
////            minAngle = mag(axy);
////        }
////        if (mag(axz) >= minAngle)
////        {
////            tmpA.xx() = mag(a.zx()); tmpA.xy() = mag(a.zy()); tmpA.xz() = mag(a.zz());
////        }
////
////        minAngle = 0;
////
////        if (mag(ayx) >= minAngle)
////        {
////            tmpA.yx() = mag(a.xx()); tmpA.yy() = mag(a.xy()); tmpA.yz() = mag(a.xz());
////            minAngle = mag(ayx);
////        }
////        if (mag(ayy) >= minAngle)
////        {
////            tmpA.yx() = mag(a.yx()); tmpA.yy() = mag(a.yy()); tmpA.yz() = mag(a.yz());
////            minAngle = mag(ayy);
////        }
////        if (mag(ayz) >= minAngle)
////        {
////            tmpA.yx() = mag(a.zx()); tmpA.yy() = mag(a.zy()); tmpA.yz() = mag(a.zz());
////        }
////
////        minAngle = 0;
////
////        if (mag(azx) >= minAngle)
////        {
////            tmpA.zx() = mag(a.xx()); tmpA.zy() = mag(a.xy()); tmpA.zz() = mag(a.xz());
////            minAngle = mag(azx);
////        }
////        if (mag(azy) >= minAngle)
////        {
////            tmpA.zx() = mag(a.yx()); tmpA.zy() = mag(a.yy()); tmpA.zz() = mag(a.yz());
////            minAngle = mag(azy);
////        }
////        if (mag(azz) >= minAngle)
////        {
////            tmpA.zx() = mag(a.zx()); tmpA.zy() = mag(a.zy()); tmpA.zz() = mag(a.zz());
////        }
////
////        alignments[aI] = tmpA;
////    }
//
//    scalar nearest = 0;
//
////    Info<< nl << "Point " << pt << endl;
////
////    FixedList<quaternion, 4> qAlignments;
////    forAll(qAlignments, aI)
////    {
//////        Info<< "    Direction " << aI << endl;
//////        Info<< "        Rot tensor" << alignments[aI] << endl;
////        qAlignments[aI] = quaternion(alignments[aI]);
////        qAlignments[aI].normalize();
//////        Info<< "        Quaternion: " << qAlignments[aI] << endl;
//////        Info<< "        Rot tensor from quat: " << qAlignments[aI].R() << endl;
////    }
//
//    tensor alignment = Foam::tensor::zero;
//    forAll(bary, pI)
//    {
////        alignment += bary[pI]*ch->vertex(pI)->alignment();
//
////        alignment += bary[pI]*alignments[pI];
//
//        // Try slerp with quaternions
//
//        // Find nearest point
//        if (bary[pI] > nearest)
//        {
//            alignment = ch->vertex(pI)->alignment();
//            nearest = bary[pI];
//        }
//    }
//
////    quaternion alignment;
//
////    alignment = qAlignments[0]*bary[0]
////              + qAlignments[1]*bary[1]
////              + qAlignments[2]*bary[2]
////              + qAlignments[3]*bary[3];
//
////    alignment = slerp(qAlignments[0], qAlignments[1], bary[0]+bary[1]+bary[2]);
////    alignment = slerp(alignment, qAlignments[2], bary[0]+bary[1]+bary[2]);
////    alignment = slerp(alignment, qAlignments[3], bary[0]+bary[1]+bary[2]);
////    alignment = slerp(alignment, qAlignments[0], bary[0]/(bary[0]+bary[1]+bary[2]));
//
////    Info<< "    Interp alignment : " << alignment << endl;
////    Info<< "    Interp rot tensor: " << alignment.R() << endl;
//
//    return alignment;
//}
//
//
//void Foam::surfaceControl::cellSizeAndAlignment
//(
//    const point& pt,
//    scalar& size,
//    tensor& alignment
//) const
//{
//    scalarList bary;
//    Cell_handle ch;
//
//    triangulatedMesh_().barycentricCoords(pt, bary, ch);
//
//    size = 0;
//    forAll(bary, pI)
//    {
//        size += bary[pI]*ch->vertex(pI)->size();
//    }
//
////    alignment = Foam::tensor::zero;
////    forAll(bary, pI)
////    {
////        alignment += bary[pI]*ch->vertex(pI)->alignment();
////    }
//
//    alignment = cellAlignment(pt);
//}


void Foam::searchableSurfaceControl::initialVertices
(
    pointField& pts,
    scalarField& sizes,
    Field<triad>& alignments
) const
{
    pts = searchableSurface_.points();

    const scalar nearFeatDistSqrCoeff = 1e-8;

    sizes.setSize(pts.size());
    alignments.setSize(pts.size());

    forAll(pts, pI)
    {
        // Is the point in the extendedFeatureEdgeMesh? If so get the
        // point normal, otherwise get the surface normal from
        // searchableSurface

        pointIndexHit info;
        label infoFeature;
        allGeometry_.findFeaturePointNearest
        (
            pts[pI],
            nearFeatDistSqrCoeff,
            info,
            infoFeature
        );

        autoPtr<triad> pointAlignment;

        if (info.hit())
        {
            const extendedFeatureEdgeMesh& features =
                allGeometry_.features()[infoFeature];

            vectorField norms = features.featurePointNormals(info.index());

            // Create a triad from these norms.
            pointAlignment.set(new triad());
            forAll(norms, nI)
            {
                pointAlignment() += norms[nI];
            }

            pointAlignment().normalize();
            pointAlignment().orthogonalize();
        }
        else
        {
            allGeometry_.findEdgeNearest
            (
                pts[pI],
                nearFeatDistSqrCoeff,
                info,
                infoFeature
            );

            if (info.hit())
            {
                const extendedFeatureEdgeMesh& features =
                    allGeometry_.features()[infoFeature];

                vectorField norms = features.edgeNormals(info.index());

                // Create a triad from these norms.
                pointAlignment.set(new triad());
                forAll(norms, nI)
                {
                    pointAlignment() += norms[nI];
                }

                pointAlignment().normalize();
                pointAlignment().orthogonalize();
            }
            else
            {
                pointField ptField(1, pts[pI]);
                scalarField distField(1, nearFeatDistSqrCoeff);
                List<pointIndexHit> infoList(1, pointIndexHit());

                searchableSurface_.findNearest(ptField, distField, infoList);

                vectorField normals(1);
                searchableSurface_.getNormal(infoList, normals);

                pointAlignment.set(new triad(normals[0]));
            }
        }

        if (!cellSizeFunction_().cellSize(pts[pI], sizes[pI]))
        {
            FatalErrorIn
            (
                "Foam::searchableSurfaceControl::initialVertices"
                "(pointField&, scalarField&, tensorField&)"
            )   << "Could not calculate cell size"
                << abort(FatalError);
        }

        alignments[pI] = pointAlignment();
    }
}


// ************************************************************************* //
