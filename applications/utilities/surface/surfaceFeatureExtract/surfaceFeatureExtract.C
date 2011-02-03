/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2011 OpenCFD Ltd.
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

Application
    surfaceFeatureExtract

Description
    Extracts and writes surface features to file

\*---------------------------------------------------------------------------*/


#include "triangle.H"
#include "triSurface.H"
#include "argList.H"
#include "Time.H"
#include "surfaceFeatures.H"
#include "featureEdgeMesh.H"
#include "treeBoundBox.H"
#include "meshTools.H"
#include "OFstream.H"
#include "triSurfaceMesh.H"
#include "vtkSurfaceWriter.H"
#include "triSurfaceFields.H"
#include "unitConversion.H"
#include "indexedOctree.H"
#include "treeDataEdge.H"

#include "buildCGALPolyhedron.H"
#include "CGALPolyhedronRings.H"
#include <CGAL/Monge_via_jet_fitting.h>
#include <CGAL/Lapack/Linear_algebra_lapack.h>
#include <CGAL/property_map.h>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void dumpBox(const treeBoundBox& bb, const fileName& fName)
{
    OFstream str(fName);

    Info<< "Dumping bounding box " << bb << " as lines to obj file "
        << str.name() << endl;


    pointField boxPoints(bb.points());

    forAll(boxPoints, i)
    {
        meshTools::writeOBJ(str, boxPoints[i]);
    }

    forAll(treeBoundBox::edges, i)
    {
        const edge& e = treeBoundBox::edges[i];

        str<< "l " << e[0]+1 <<  ' ' << e[1]+1 << nl;
    }
}


// Deletes all edges inside/outside bounding box from set.
void deleteBox
(
    const triSurface& surf,
    const treeBoundBox& bb,
    const bool removeInside,
    List<surfaceFeatures::edgeStatus>& edgeStat
)
{
    forAll(edgeStat, edgeI)
    {
        const point eMid = surf.edges()[edgeI].centre(surf.localPoints());

        if (removeInside ? bb.contains(eMid) : !bb.contains(eMid))
        {
            edgeStat[edgeI] = surfaceFeatures::NONE;
        }
    }
}


void drawHitProblem
(
    label fI,
    const triSurface& surf,
    const pointField& start,
    const pointField& faceCentres,
    const pointField& end,
    const List<pointIndexHit>& hitInfo
)
{
    Info<< nl << "# findLineAll did not hit its own face."
        << nl << "# fI " << fI
        << nl << "# start " << start[fI]
        << nl << "# f centre " << faceCentres[fI]
        << nl << "# end " << end[fI]
        << nl << "# hitInfo " << hitInfo
        << endl;

    meshTools::writeOBJ(Info, start[fI]);
    meshTools::writeOBJ(Info, faceCentres[fI]);
    meshTools::writeOBJ(Info, end[fI]);

    Info<< "l 1 2 3" << endl;

    meshTools::writeOBJ(Info, surf.points()[surf[fI][0]]);
    meshTools::writeOBJ(Info, surf.points()[surf[fI][1]]);
    meshTools::writeOBJ(Info, surf.points()[surf[fI][2]]);

    Info<< "f 4 5 6" << endl;

    forAll(hitInfo, hI)
    {
        label hFI = hitInfo[hI].index();

        meshTools::writeOBJ(Info, surf.points()[surf[hFI][0]]);
        meshTools::writeOBJ(Info, surf.points()[surf[hFI][1]]);
        meshTools::writeOBJ(Info, surf.points()[surf[hFI][2]]);

        Info<< "f "
            << 3*hI + 7 << " "
            << 3*hI + 8 << " "
            << 3*hI + 9
            << endl;
    }
}


scalarField curvature(const triSurface& surf)
{
    scalarField k(surf.points().size(), 0);

    Polyhedron P;

    buildCGALPolyhedron convert(surf);
    P.delegate(convert);

    // Info<< "Created CGAL Polyhedron with " << label(P.size_of_vertices())
    //     << " vertices and " << label(P.size_of_facets())
    //     << " facets. " << endl;

    // The rest of this function adapted from
    //     CGAL-3.7/examples/Jet_fitting_3/Mesh_estimation.cpp
    // Licensed under CGAL-3.7/LICENSE.FREE_USE

    // Copyright (c) 1996,1997,1998,1999,2000,2001,2002,2003,2004,2005,2006,2007
    // Utrecht University (The Netherlands), ETH Zurich (Switzerland), Freie
    // Universitaet Berlin (Germany), INRIA Sophia-Antipolis (France),
    // Martin-Luther-University Halle-Wittenberg (Germany), Max-Planck-Institute
    // Saarbruecken (Germany), RISC Linz (Austria), and Tel-Aviv University
    // (Israel).  All rights reserved.

    // Permission is hereby granted, free of charge, to any person obtaining a
    // copy of this software and associated documentation files (the
    // "Software"), to deal in the Software without restriction, including
    // without limitation the rights to use, copy, modify, merge, publish,
    // distribute, sublicense, and/or sell copies of the Software, and to permit
    // persons to whom the Software is furnished to do so, subject to the
    // following conditions:

    // THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
    // OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
    // MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
    // IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
    // CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT
    // OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
    // THE USE OR OTHER DEALINGS IN THE SOFTWARE.

     //Vertex property map, with std::map
    typedef std::map<Vertex*, int> Vertex2int_map_type;
    typedef boost::associative_property_map< Vertex2int_map_type >
        Vertex_PM_type;
    typedef T_PolyhedralSurf_rings<Polyhedron, Vertex_PM_type > Poly_rings;

    typedef CGAL::Monge_via_jet_fitting<Kernel>         Monge_via_jet_fitting;
    typedef Monge_via_jet_fitting::Monge_form           Monge_form;

    std::vector<Point_3> in_points;  //container for data points

    // default parameter values and global variables
    unsigned int d_fitting = 2;
    unsigned int d_monge = 2;
    unsigned int min_nb_points = (d_fitting + 1)*(d_fitting + 2)/2;

    //initialize the tag of all vertices to -1
    Vertex_iterator vitb = P.vertices_begin();
    Vertex_iterator vite = P.vertices_end();

    Vertex2int_map_type vertex2props;
    Vertex_PM_type vpm(vertex2props);

    CGAL_For_all(vitb, vite)
    {
        put(vpm, &(*vitb), -1);
    }

    vite = P.vertices_end();

    label vertI = 0;

    for (vitb = P.vertices_begin(); vitb != vite; vitb++)
    {
        //initialize
        Vertex* v = &(*vitb);

        //gather points around the vertex using rings
        // From: gather_fitting_points(v, in_points, vpm);
        {
            std::vector<Vertex*> gathered;
            in_points.clear();

            Poly_rings::collect_enough_rings(v, min_nb_points, gathered, vpm);

            //store the gathered points
            std::vector<Vertex*>::iterator itb = gathered.begin();
            std::vector<Vertex*>::iterator ite = gathered.end();

            CGAL_For_all(itb, ite)
            {
                in_points.push_back((*itb)->point());
            }
        }

        //skip if the nb of points is to small
        if ( in_points.size() < min_nb_points )
        {
            std::cerr
                << "not enough pts for fitting this vertex"
                << in_points.size()
                << std::endl;

            continue;
        }

        // perform the fitting
        Monge_via_jet_fitting monge_fit;

        Monge_form monge_form = monge_fit
        (
            in_points.begin(),
            in_points.end(),
            d_fitting,
            d_monge
        );

        // std::cout<< monge_form << std::endl;

        // std::cout<< "k1 " << monge_form.principal_curvatures(0) << std::endl;
        // std::cout<< "k2 " << monge_form.principal_curvatures(1) << std::endl;

        // std::vector<Point_3>::iterator itbp = in_points.begin();
        // std::vector<Point_3>::iterator itep = in_points.end();

        // std::cout << "in_points list : " << std::endl;

        // for (; itbp != itep; itbp++)
        // {
        //     std::cout << *itbp << std::endl;
        // }

        // std::cout << "--- vertex " <<  vertI
        //     << " : " << v->point() << std::endl
        //     << "number of points used : " << in_points.size()
        //     << std::endl;

        k[vertI++] = Foam::sqrt
        (
            sqr(monge_form.principal_curvatures(0))
          + sqr(monge_form.principal_curvatures(1))
        );
    }

    return k;
}


// Main program:

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "extract and write surface features to file"
    );
    argList::noParallel();
    argList::validArgs.append("surface");
    argList::validArgs.append("output set");

    argList::addOption
    (
        "includedAngle",
        "degrees",
        "construct feature set from included angle [0..180]"
    );
    argList::addOption
    (
        "set",
        "name",
        "use existing feature set from file"
    );
    argList::addOption
    (
        "minLen",
        "scalar",
        "remove features shorter than the specified cumulative length"
    );
    argList::addOption
    (
        "minElem",
        "int",
        "remove features with fewer than the specified number of edges"
    );
    argList::addOption
    (
        "subsetBox",
        "((x0 y0 z0)(x1 y1 z1))",
        "remove edges outside specified bounding box"
    );
    argList::addOption
    (
        "deleteBox",
        "((x0 y0 z0)(x1 y1 z1))",
        "remove edges within specified bounding box"
    );
    argList::addBoolOption
    (
        "writeObj",
        "write featureEdgeMesh obj files"
    );
    argList::addOption
    (
        "closeness",
        "scalar",
        "span to look for surface closeness"
    );
    argList::addOption
    (
        "featureProximity",
        "scalar",
        "distance to look for close features"
    );
    argList::addBoolOption
    (
        "writeVTK",
        "write surface property VTK files"
    );

#   include "setRootCase.H"
#   include "createTime.H"

    Info<< "Feature line extraction is only valid on closed manifold surfaces."
        << endl;

    bool writeVTK = args.optionFound("writeVTK");
    bool writeObj = args.optionFound("writeObj");

    const fileName surfFileName = args[1];
    const fileName outFileName  = args[2];

    Info<< "Surface            : " << surfFileName << nl
        << "Output feature set : " << outFileName << nl
        << endl;

    fileName sFeatFileName = surfFileName.lessExt().name();


    // Read
    // ~~~~

    triSurface surf(surfFileName);

    Info<< "Statistics:" << endl;
    surf.writeStats(Info);
    Info<< endl;

    faceList faces(surf.size());

    forAll(surf, fI)
    {
        faces[fI] = surf[fI].triFaceFace();
    }


    // Either construct features from surface&featureangle or read set.
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    surfaceFeatures set(surf);

    if (args.optionFound("set"))
    {
        const fileName setName = args["set"];

        Info<< "Reading existing feature set from file " << setName << endl;

        set = surfaceFeatures(surf, setName);
    }
    else if (args.optionFound("includedAngle"))
    {
        const scalar includedAngle = args.optionRead<scalar>("includedAngle");

        Info<< "Constructing feature set from included angle " << includedAngle
            << endl;

        set = surfaceFeatures(surf, includedAngle);

        // Info<< nl << "Writing initial features" << endl;
        // set.write("initial.fSet");
        // set.writeObj("initial");
    }
    else
    {
        FatalErrorIn(args.executable())
            << "No initial feature set. Provide either one"
            << " of -set (to read existing set)" << nl
            << " or -includedAngle (to new set construct from angle)"
            << exit(FatalError);
    }

    Info<< nl
        << "Initial feature set:" << nl
        << "    feature points : " << set.featurePoints().size() << nl
        << "    feature edges  : " << set.featureEdges().size() << nl
        << "    of which" << nl
        << "        region edges   : " << set.nRegionEdges() << nl
        << "        external edges : " << set.nExternalEdges() << nl
        << "        internal edges : " << set.nInternalEdges() << nl
        << endl;


    // Trim set
    // ~~~~~~~~

    scalar minLen = -GREAT;
    if (args.optionReadIfPresent("minLen", minLen))
    {
        Info<< "Removing features of length < " << minLen << endl;
    }

    label minElem = 0;
    if (args.optionReadIfPresent("minElem", minElem))
    {
        Info<< "Removing features with number of edges < " << minElem << endl;
    }

    // Trim away small groups of features
    if (minElem > 0 || minLen > 0)
    {
        set.trimFeatures(minLen, minElem);
        Info<< endl << "Removed small features" << endl;
    }


    // Subset
    // ~~~~~~

    // Convert to marked edges, points
    List<surfaceFeatures::edgeStatus> edgeStat(set.toStatus());

    if (args.optionFound("subsetBox"))
    {
        treeBoundBox bb
        (
            args.optionLookup("subsetBox")()
        );

        Info<< "Removing all edges outside bb " << bb << endl;
        dumpBox(bb, "subsetBox.obj");

        deleteBox
        (
            surf,
            bb,
            false,
            edgeStat
        );
    }
    else if (args.optionFound("deleteBox"))
    {
        treeBoundBox bb
        (
            args.optionLookup("deleteBox")()
        );

        Info<< "Removing all edges inside bb " << bb << endl;
        dumpBox(bb, "deleteBox.obj");

        deleteBox
        (
            surf,
            bb,
            true,
            edgeStat
        );
    }

    surfaceFeatures newSet(surf);
    newSet.setFromStatus(edgeStat);

    Info<< endl << "Writing trimmed features to "
        << runTime.constant()/"featureEdgeMesh"/outFileName << endl;
    newSet.write(runTime.constant()/"featureEdgeMesh"/outFileName);

    // Info<< endl << "Writing edge objs." << endl;
    // newSet.writeObj("final");

    Info<< nl
        << "Final feature set:" << nl
        << "    feature points : " << newSet.featurePoints().size() << nl
        << "    feature edges  : " << newSet.featureEdges().size() << nl
        << "    of which" << nl
        << "        region edges   : " << newSet.nRegionEdges() << nl
        << "        external edges : " << newSet.nExternalEdges() << nl
        << "        internal edges : " << newSet.nInternalEdges() << nl
        << endl;

    // Extracting and writing a featureEdgeMesh

    Pout<< nl << "Writing featureEdgeMesh to constant/featureEdgeMesh."
        << endl;

    featureEdgeMesh feMesh
    (
        IOobject
        (
            sFeatFileName + ".featureEdgeMesh",
            runTime.constant(),
            "featureEdgeMesh",
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        newSet
    );

    feMesh.write();

    if (writeObj)
    {
        feMesh.writeObj(runTime.constant()/"featureEdgeMesh"/sFeatFileName);
    };

    triSurfaceMesh searchSurf
    (
        IOobject
        (
            sFeatFileName + ".closeness",
            runTime.constant(),
            "featureEdgeMesh",
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        surf
    );

    // Find close features

    // // Dummy trim operation to mark features
    // labelList featureEdgeIndexing = newSet.trimFeatures(-GREAT, 0);

    // scalarField surfacePtFeatureIndex(surf.points().size(), -1);

    // forAll(newSet.featureEdges(), eI)
    // {
    //     const edge& e = surf.edges()[newSet.featureEdges()[eI]];

    //     surfacePtFeatureIndex[surf.meshPoints()[e.start()]] =
    //     featureEdgeIndexing[newSet.featureEdges()[eI]];

    //     surfacePtFeatureIndex[surf.meshPoints()[e.end()]] =
    //     featureEdgeIndexing[newSet.featureEdges()[eI]];
    // }

    // if (writeVTK)
    // {
    //     vtkSurfaceWriter().write
    //     (
    //         runTime.constant()/"triSurface",    // outputDir
    //         sFeatFileName,                      // surfaceName
    //         surf.points(),
    //         faces,
    //         "surfacePtFeatureIndex",            // fieldName
    //         surfacePtFeatureIndex,
    //         true,                               // isNodeValues
    //         true                                // verbose
    //     );
    // }

    // Random rndGen(343267);

    // treeBoundBox surfBB
    // (
    //     treeBoundBox(searchSurf.bounds()).extend(rndGen, 1e-4)
    // );

    // surfBB.min() -= Foam::point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
    // surfBB.max() += Foam::point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

    // indexedOctree<treeDataEdge> ftEdTree
    // (
    //     treeDataEdge
    //     (
    //         false,
    //         surf.edges(),
    //         surf.localPoints(),
    //         newSet.featureEdges()
    //     ),
    //     surfBB,
    //     8,      // maxLevel
    //     10,     // leafsize
    //     3.0     // duplicity
    // );

    // labelList nearPoints = ftEdTree.findBox
    // (
    //     treeBoundBox
    //     (
    //         sPt - featureSearchSpan*Foam::vector::one,
    //         sPt + featureSearchSpan*Foam::vector::one
    //     )
    // );

    // Examine curvature, feature proximity and internal and external closeness.

    // Internal and external closeness

    // Prepare start and end points for intersection tests

    const vectorField& normals = searchSurf.faceNormals();

    scalar span = searchSurf.bounds().mag();

    args.optionReadIfPresent("closeness", span);

    scalar externalAngleTolerance = 10;
    scalar externalToleranceCosAngle = Foam::cos
    (
        degToRad(180 - externalAngleTolerance)
    );

    scalar internalAngleTolerance = 45;
    scalar internalToleranceCosAngle = Foam::cos
    (
        degToRad(180 - internalAngleTolerance)
    );

    Info<< "externalToleranceCosAngle: " << externalToleranceCosAngle << nl
        << "internalToleranceCosAngle: " << internalToleranceCosAngle
        << endl;

    // Info<< "span " << span << endl;

    pointField start = searchSurf.faceCentres() - span*normals;
    pointField end = searchSurf.faceCentres() + span*normals;
    const pointField& faceCentres = searchSurf.faceCentres();

    List<List<pointIndexHit> > allHitInfo;

    // Find all intersections (in order)
    searchSurf.findLineAll(start, end, allHitInfo);

    scalarField internalCloseness(start.size(), GREAT);
    scalarField externalCloseness(start.size(), GREAT);

    forAll(allHitInfo, fI)
    {
        const List<pointIndexHit>& hitInfo = allHitInfo[fI];

        if (hitInfo.size() < 1)
        {
            drawHitProblem(fI, surf, start, faceCentres, end, hitInfo);

            // FatalErrorIn(args.executable())
            //     << "findLineAll did not hit its own face."
            //     << exit(FatalError);
        }
        else if (hitInfo.size() == 1)
        {
            if (!hitInfo[0].hit())
            {
                // FatalErrorIn(args.executable())
                //     << "findLineAll did not hit any face."
                //     << exit(FatalError);
            }
            else if (hitInfo[0].index() != fI)
            {
                drawHitProblem(fI, surf, start, faceCentres, end, hitInfo);

                // FatalErrorIn(args.executable())
                //     << "findLineAll did not hit its own face."
                //     << exit(FatalError);
            }
        }
        else
        {
            label ownHitI = -1;

            forAll(hitInfo, hI)
            {
                // Find the hit on the triangle that launched the ray

                if (hitInfo[hI].index() == fI)
                {
                    ownHitI = hI;

                    break;
                }
            }

            if (ownHitI < 0)
            {
                drawHitProblem(fI, surf, start, faceCentres, end, hitInfo);

                // FatalErrorIn(args.executable())
                //     << "findLineAll did not hit its own face."
                //     << exit(FatalError);
            }
            else if (ownHitI == 0)
            {
                // There are no internal hits, the first hit is the closest
                // external hit

                if
                (
                    (normals[fI] & normals[hitInfo[ownHitI + 1].index()])
                  < externalToleranceCosAngle
                )
                {
                    externalCloseness[fI] = mag
                    (
                        faceCentres[fI] - hitInfo[ownHitI + 1].hitPoint()
                    );
                }
            }
            else if (ownHitI == hitInfo.size() - 1)
            {
                // There are no external hits, the last but one hit is the
                // closest internal hit

                if
                (
                    (normals[fI] & normals[hitInfo[ownHitI - 1].index()])
                  < internalToleranceCosAngle
                )
                {
                    internalCloseness[fI] = mag
                    (
                        faceCentres[fI] - hitInfo[ownHitI - 1].hitPoint()
                    );
                }
            }
            else
            {
                if
                (
                    (normals[fI] & normals[hitInfo[ownHitI + 1].index()])
                  < externalToleranceCosAngle
                )
                {
                    externalCloseness[fI] = mag
                    (
                        faceCentres[fI] - hitInfo[ownHitI + 1].hitPoint()
                    );
                }

                if
                (
                    (normals[fI] & normals[hitInfo[ownHitI - 1].index()])
                  < internalToleranceCosAngle
                )
                {
                    internalCloseness[fI] = mag
                    (
                        faceCentres[fI] - hitInfo[ownHitI - 1].hitPoint()
                    );
                }
            }
        }
    }

    triSurfaceScalarField internalClosenessField
    (
        IOobject
        (
            sFeatFileName + ".internalCloseness",
            runTime.constant(),
            "featureEdgeMesh",
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        surf,
        dimLength,
        internalCloseness
    );

    internalClosenessField.write();

    triSurfaceScalarField externalClosenessField
    (
        IOobject
        (
            sFeatFileName + ".externalCloseness",
            runTime.constant(),
            "featureEdgeMesh",
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        surf,
        dimLength,
        externalCloseness
    );

    externalClosenessField.write();

    scalarField k = curvature(surf);

    // Modify the curvature values on feature edges and points to be zero.

    forAll(newSet.featureEdges(), fEI)
    {
        const edge& e = surf.edges()[newSet.featureEdges()[fEI]];

        k[surf.meshPoints()[e.start()]] = 0.0;
        k[surf.meshPoints()[e.end()]] = 0.0;
    }

    triSurfacePointScalarField kField
    (
        IOobject
        (
            sFeatFileName + ".curvature",
            runTime.constant(),
            "featureEdgeMesh",
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        surf,
        dimLength,
        k
    );

    kField.write();

    if (writeVTK)
    {
        vtkSurfaceWriter().write
        (
            runTime.constant()/"triSurface",    // outputDir
            sFeatFileName,                      // surfaceName
            surf.points(),
            faces,
            "internalCloseness",                // fieldName
            internalCloseness,
            false,                              // isNodeValues
            true                                // verbose
        );

        vtkSurfaceWriter().write
        (
            runTime.constant()/"triSurface",    // outputDir
            sFeatFileName,                      // surfaceName
            surf.points(),
            faces,
            "externalCloseness",                // fieldName
            externalCloseness,
            false,                              // isNodeValues
            true                                // verbose
        );

        vtkSurfaceWriter().write
        (
            runTime.constant()/"triSurface",    // outputDir
            sFeatFileName,                      // surfaceName
            surf.points(),
            faces,
            "curvature",                        // fieldName
            k,
            true,                               // isNodeValues
            true                                // verbose
        );
    }

    Info<< nl << "End" << endl;

    return 0;
}


// ************************************************************************* //
