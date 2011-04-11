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
#include "extendedFeatureEdgeMesh.H"
#include "treeBoundBox.H"
#include "meshTools.H"
#include "OFstream.H"
#include "unitConversion.H"

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


// Unmark non-manifold edges if individual triangles are not features
void unmarkBaffles
(
    const triSurface& surf,
    const scalar includedAngle,
    List<surfaceFeatures::edgeStatus>& edgeStat
)
{
    scalar minCos = Foam::cos(degToRad(180.0 - includedAngle));

    const labelListList& edgeFaces = surf.edgeFaces();

    forAll(edgeFaces, edgeI)
    {
        const labelList& eFaces = edgeFaces[edgeI];

        if (eFaces.size() > 2)
        {
            label i0 = eFaces[0];
            //const labelledTri& f0 = surf[i0];
            const vector& n0 = surf.faceNormals()[i0];

            //Pout<< "edge:" << edgeI << " n0:" << n0 << endl;

            bool same = true;

            for (label i = 1; i < eFaces.size(); i++)
            {
                //const labelledTri& f = surf[i];
                const vector& n = surf.faceNormals()[eFaces[i]];

                //Pout<< "    mag(n&n0): " << mag(n&n0) << endl;

                if (mag(n&n0) < minCos)
                {
                    same = false;
                    break;
                }
            }

            if (same)
            {
                edgeStat[edgeI] = surfaceFeatures::NONE;
            }
        }
    }
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
        "write extendedFeatureEdgeMesh obj files"
    );

#   include "setRootCase.H"
#   include "createTime.H"

    Info<< "Feature line extraction is only valid on closed manifold surfaces."
        << endl;

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

    Info<< endl << "Writing trimmed features to " << outFileName << endl;
    newSet.write(outFileName);

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

    // Extracting and writing a extendedFeatureEdgeMesh

    extendedFeatureEdgeMesh feMesh
    (
        newSet,
        runTime,
        sFeatFileName + ".extendedFeatureEdgeMesh"
    );

    Info<< nl << "Writing extendedFeatureEdgeMesh to " << feMesh.objectPath()
        << endl;

    if (writeObj)
    {
        feMesh.writeObj(surfFileName.lessExt().name());
    }

    feMesh.write();


    // Write a featureEdgeMesh for backwards compatibility
    {
        featureEdgeMesh bfeMesh
        (
            IOobject
            (
                surfFileName.lessExt().name() + ".eMesh",   // name
                runTime.constant(), // instance
                "triSurface",
                runTime,            // registry
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

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
