/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2017 OpenCFD Ltd.
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

Group
    grpSurfaceUtilities

Description
    Extracts and writes surface features to file. All but the basic feature
    extraction is a work-in-progress.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "triSurface.H"
#include "triSurfaceTools.H"
#include "edgeMeshTools.H"
#include "surfaceFeaturesExtraction.H"
#include "surfaceIntersection.H"
#include "featureEdgeMesh.H"
#include "extendedFeatureEdgeMesh.H"
#include "treeBoundBox.H"
#include "meshTools.H"
#include "OBJstream.H"
#include "triSurfaceMesh.H"
#include "vtkSurfaceWriter.H"
#include "unitConversion.H"
#include "plane.H"
#include "point.H"
#include "triSurfaceLoader.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "extract and write surface features to file"
    );
    argList::noParallel();
    argList::noFunctionObjects();

    #include "addDictOption.H"

    #include "setRootCase.H"
    #include "createTime.H"

    const word dictName("surfaceFeatureExtractDict");
    #include "setSystemRunTimeDictionaryIO.H"

    Info<< "Reading " << dictName << nl << endl;
    const IOdictionary dict(dictIO);

    // Loader for available triSurface surface files
    triSurfaceLoader loader(runTime);

    // Where to write VTK output files
    const fileName vtkOutputDir = runTime.constantPath()/"triSurface";

    forAllConstIter(dictionary, dict, iter)
    {
        const word& dictName = iter().keyword();

        if (!iter().isDict())
        {
            continue;
        }
        const dictionary& surfaceDict = iter().dict();

        if (!surfaceDict.found("extractionMethod"))
        {
            continue;
        }

        autoPtr<surfaceFeaturesExtraction::method> extractor =
            surfaceFeaturesExtraction::method::New
            (
                surfaceDict
            );

        // The output name, cleansed of extensions
        // Optional "output" entry, or the dictionary name.
        const word outputName =
            fileName
            (
                surfaceDict.lookupOrDefault<word>("output", dictName)
            ).lessExt();

        // The "surfaces" entry is normally optional, but if the sub-dictionary
        // is itself called "surfaces", then this becomes mandatory.
        // This provides a simple means of handling both situations without an
        // additional switch.
        if
        (
            dictName == "surfaces"  // mandatory
         || surfaceDict.found("surfaces")   // or optional
        )
        {
            loader.select(wordReList(surfaceDict.lookup("surfaces")));
        }
        else
        {
            loader.select(dictName);
        }

        if (loader.selected().empty())
        {
            FatalErrorInFunction
                << "No surfaces specified/found for entry: "
                << dictName << exit(FatalError);
        }
        // DebugVar(loader.available());
        // DebugVar(outputName);


        Info<< "Surfaces           : ";
        if (loader.selected().size() == 1)
        {
            Info<< loader.selected()[0] << nl;
        }
        else
        {
            Info<< flatOutput(loader.selected()) << nl;
        }
        Info<< "Output             : " << outputName << nl;

        triSurfaceLoader::loadingOption loadingOption =
            triSurfaceLoader::loadingOptionNames.lookupOrDefault
            (
                "loadingOption",
                surfaceDict,
                triSurfaceLoader::loadingOption::OFFSET_REGION
            );

        Info<<"loading with "
            << triSurfaceLoader::loadingOptionNames[loadingOption]
            << endl;


        // Load a single file, or load and combine multiple selected files
        autoPtr<triSurface> surfPtr = loader.load(loadingOption);
        if (!surfPtr.valid() || surfPtr().empty())
        {
            FatalErrorInFunction
                << "Problem loading surface(s) for entry: "
                << dictName << exit(FatalError);
        }

        triSurface surf = surfPtr();

        const Switch writeVTK = surfaceDict.lookupOrDefault<Switch>
        (
            "writeVTK",
            Switch::OFF
        );
        const Switch writeObj = surfaceDict.lookupOrDefault<Switch>
        (
            "writeObj",
            Switch::OFF
        );

        Info<< "write VTK: " <<  writeVTK << nl;

        Info<< "Feature line extraction is only valid on closed manifold "
            << "surfaces." << nl;

        Info<< nl << "Statistics:" << nl;
        surf.writeStats(Info);
        Info<< nl;

        // Need a copy as plain faces if outputting VTK format
        faceList faces;
        if (writeVTK)
        {
            faces.setSize(surf.size());
            forAll(surf, fi)
            {
                faces[fi] = surf[fi].triFaceFace();
            }
        }

        //
        // Extract features using the preferred extraction method
        //
        autoPtr<surfaceFeatures> features = extractor().features(surf);

        // Trim set
        // ~~~~~~~~

        // Option: "trimFeatures" (dictionary)
        if (surfaceDict.isDict("trimFeatures"))
        {
            const dictionary& trimDict = surfaceDict.subDict("trimFeatures");

            const scalar minLen =
                trimDict.lookupOrDefault<scalar>("minLen", 0);
            const label minElem =
                trimDict.lookupOrDefault<label>("minElem", 0);

            // Trim away small groups of features
            if (minLen > 0 || minElem > 0)
            {
                if (minLen > 0)
                {
                    Info<< "Removing features of length < "
                        << minLen << endl;
                }
                if (minElem > 0)
                {
                    Info<< "Removing features with number of edges < "
                        << minElem << endl;
                }

                features().trimFeatures
                (
                    minLen, minElem, extractor().includedAngle()
                );
            }
        }

        // Subset
        // ~~~~~~

        // Convert to marked edges, points
        List<surfaceFeatures::edgeStatus> edgeStat(features().toStatus());

        // Option: "subsetFeatures" (dictionary)
        if (surfaceDict.isDict("subsetFeatures"))
        {
            const dictionary& subsetDict = surfaceDict.subDict
            (
                "subsetFeatures"
            );

            // Suboption: "insideBox"
            if (subsetDict.found("insideBox"))
            {
                treeBoundBox bb(subsetDict.lookup("insideBox")());

                Info<< "Subset edges inside box " << bb << endl;
                features().subsetBox(edgeStat, bb);

                {
                    OBJstream os("subsetBox.obj");

                    Info<< "Dumping bounding box " << bb
                        << " as lines to obj file "
                        << os.name() << endl;

                    os.write(bb);
                }
            }
            // Suboption: "outsideBox"
            else if (subsetDict.found("outsideBox"))
            {
                treeBoundBox bb(subsetDict.lookup("outsideBox")());

                Info<< "Exclude edges outside box " << bb << endl;
                features().excludeBox(edgeStat, bb);

                {
                    OBJstream os("deleteBox.obj");

                    Info<< "Dumping bounding box " << bb
                        << " as lines to obj file "
                        << os.name() << endl;

                    os.write(bb);
                }
            }

            // Suboption: "nonManifoldEdges" (false: remove non-manifold edges)
            if (!subsetDict.lookupOrDefault<bool>("nonManifoldEdges", true))
            {
                Info<< "Removing all non-manifold edges"
                    << " (edges with > 2 connected faces) unless they"
                    << " cross multiple regions" << endl;

                features().checkFlatRegionEdge
                (
                    edgeStat,
                    1e-5,   // tol
                    extractor().includedAngle()
                );
            }

            // Suboption: "openEdges" (false: remove open edges)
            if (!subsetDict.lookupOrDefault<bool>("openEdges", true))
            {
                Info<< "Removing all open edges"
                    << " (edges with 1 connected face)" << endl;

                features().excludeOpen(edgeStat);
            }

            // Suboption: "plane"
            if (subsetDict.found("plane"))
            {
                plane cutPlane(subsetDict.lookup("plane")());

                Info<< "Only include feature edges that intersect the plane"
                    << " with normal " << cutPlane.normal()
                    << " and base point " << cutPlane.refPoint() << endl;

                features().subsetPlane(edgeStat, cutPlane);
            }
        }

        surfaceFeatures newSet(surf);
        newSet.setFromStatus(edgeStat, extractor().includedAngle());

        Info<< nl << "Initial ";
        newSet.writeStats(Info);

        boolList surfBaffleRegions(surf.patches().size(), false);
        if (surfaceDict.found("baffles"))
        {
            wordReList baffleSelect(surfaceDict.lookup("baffles"));

            wordList patchNames(surf.patches().size());
            forAll(surf.patches(), patchi)
            {
                patchNames[patchi] = surf.patches()[patchi].name();
            }

            labelList indices = findStrings(baffleSelect, patchNames);
            forAll(indices, patchi)
            {
                surfBaffleRegions[patchi] = true;
            }

            if (indices.size())
            {
                Info<< "Adding " << indices.size() << " baffle regions: (";

                forAll(surfBaffleRegions, patchi)
                {
                    if (surfBaffleRegions[patchi])
                    {
                        Info<< ' ' << patchNames[patchi];
                    }
                }
                Info<< " )" << nl << nl;
            }
        }

        // Extracting and writing a extendedFeatureEdgeMesh
        extendedFeatureEdgeMesh feMesh
        (
            newSet,
            runTime,
            outputName + ".extendedFeatureEdgeMesh",
            surfBaffleRegions
        );


        if (surfaceDict.isDict("addFeatures"))
        {
            const word addFeName = surfaceDict.subDict("addFeatures")["name"];
            Info<< "Adding (without merging) features from " << addFeName
                << nl << endl;

            extendedFeatureEdgeMesh addFeMesh
            (
                IOobject
                (
                    addFeName,
                    runTime.time().constant(),
                    "extendedFeatureEdgeMesh",
                    runTime.time(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            );
            Info<< "Read " << addFeMesh.name() << nl;
            edgeMeshTools::writeStats(Info, addFeMesh);

            feMesh.add(addFeMesh);
        }

        const surfaceIntersection::intersectionType selfIntersect =
            surfaceIntersection::selfIntersectionNames.lookupOrDefault
            (
                "intersectionMethod",
                surfaceDict,
                surfaceIntersection::NONE
            );

        if (selfIntersect != surfaceIntersection::NONE)
        {
            triSurfaceSearch query(surf);
            surfaceIntersection intersect(query, surfaceDict);

            // Remove rounding noise - could make adjustable
            intersect.mergePoints(10*SMALL);

            labelPair sizeInfo
            (
                intersect.cutPoints().size(),
                intersect.cutEdges().size()
            );

            if (intersect.cutEdges().size())
            {
                extendedEdgeMesh addMesh
                (
                    intersect.cutPoints(),
                    intersect.cutEdges()
                );

                feMesh.add(addMesh);

                sizeInfo[0] = addMesh.points().size();
                sizeInfo[1] = addMesh.edges().size();
            }
            Info<< nl
                << "intersection: "
                << surfaceIntersection::selfIntersectionNames[selfIntersect]
                << nl
                << "    points : " << sizeInfo[0] << nl
                << "    edges  : " << sizeInfo[1] << nl;
        }

        Info<< nl << "Final ";
        edgeMeshTools::writeStats(Info, feMesh);

        Info<< nl << "Writing extendedFeatureEdgeMesh to "
            << feMesh.objectPath() << endl;

        mkDir(feMesh.path());

        if (writeObj)
        {
            feMesh.writeObj(feMesh.path()/outputName);
        }

        feMesh.write();

        // Write a featureEdgeMesh for backwards compatibility
        if (true)
        {
            featureEdgeMesh bfeMesh
            (
                IOobject
                (
                    outputName + ".eMesh",      // name
                    runTime.constant(),         // instance
                    "triSurface",
                    runTime,                    // registry
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE,
                    false
                ),
                feMesh.points(),
                feMesh.edges()
            );

            Info<< nl << "Writing featureEdgeMesh to "
                << bfeMesh.objectPath() << endl;

            bfeMesh.regIOobject::write();
        }

        // Option: "closeness"
        if (surfaceDict.lookupOrDefault<bool>("closeness", false))
        {
            Pair<tmp<scalarField>> tcloseness =
                triSurfaceTools::writeCloseness
                (
                    runTime,
                    outputName,
                    surf,
                    45,  // internalAngleTolerance
                    10   // externalAngleTolerance
                );

            if (writeVTK)
            {
                vtkSurfaceWriter().write
                (
                    vtkOutputDir,
                    outputName,
                    meshedSurfRef
                    (
                        surf.points(),
                        faces
                    ),
                    "internalCloseness",                // fieldName
                    tcloseness[0](),
                    false,                              // isNodeValues
                    true                                // verbose
                );

                vtkSurfaceWriter().write
                (
                    vtkOutputDir,
                    outputName,
                    meshedSurfRef
                    (
                        surf.points(),
                        faces
                    ),
                    "externalCloseness",                // fieldName
                    tcloseness[1](),
                    false,                              // isNodeValues
                    true                                // verbose
                );
            }
        }

        // Option: "curvature"
        if (surfaceDict.lookupOrDefault<bool>("curvature", false))
        {
            tmp<scalarField> tcurvatureField =
                triSurfaceTools::writeCurvature
                (
                    runTime,
                    outputName,
                    surf
                );

            if (writeVTK)
            {
                vtkSurfaceWriter().write
                (
                    vtkOutputDir,
                    outputName,
                    meshedSurfRef
                    (
                        surf.points(),
                        faces
                    ),
                    "curvature",                        // fieldName
                    tcurvatureField(),
                    true,                               // isNodeValues
                    true                                // verbose
                );
            }
        }

        // Option: "featureProximity"
        if (surfaceDict.lookupOrDefault<bool>("featureProximity", false))
        {
            tmp<scalarField> tproximity =
                edgeMeshTools::writeFeatureProximity
                (
                    runTime,
                    outputName,
                    feMesh,
                    surf,
                    readScalar(surfaceDict.lookup("maxFeatureProximity"))
                );

            if (writeVTK)
            {
                vtkSurfaceWriter().write
                (
                    vtkOutputDir,
                    outputName,
                    meshedSurfRef
                    (
                        surf.points(),
                        faces
                    ),
                    "featureProximity",                 // fieldName
                    tproximity(),
                    false,                              // isNodeValues
                    true                                // verbose
                );
            }
        }

        Info<< endl;
    }

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
