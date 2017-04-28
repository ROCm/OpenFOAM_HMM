/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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
    Test-surfaceIntersection

Description
    Test surface-surface intersection

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "triSurface.H"
#include "triSurfaceMesh.H"
#include "surfaceIntersection.H"
#include "OFstream.H"

using namespace Foam;

autoPtr<triSurface> loadSurface
(
    const Foam::Time& runTime,
    const fileName& surfName
)
{
    Info<< "Reading surface " << surfName << endl;

    const fileName fallback =
        runTime.constantPath()/triSurfaceMesh::meshSubDir/surfName;

    autoPtr<triSurface> surfPtr;
    if (isFile(surfName))
    {
        surfPtr.set(new triSurface(surfName));
    }
    else if (isFile(fallback))
    {
        surfPtr.set(new triSurface(fallback));
    }
    else
    {
        FatalErrorInFunction
            << "No such file:" << surfName << exit(FatalError);
    }

    return surfPtr;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Intersection of two surfaces. Writes obj file"
    );
    argList::addBoolOption
    (
        "debug2",
        "set surfaceIntersection debug=2"
    );
    argList::addBoolOption
    (
        "debug4",
        "set surfaceIntersection debug=4"
    );
    argList::addBoolOption
    (
        "print",
        "print information about cuts, etc"
    );
    argList::addBoolOption
    (
        "mergeEdges",
        "merge duplicate edges"
    );
    argList::addOption
    (
        "mergePoints",
        "mergeTol",
        "merge points (and edges) using the specified tolerance"
    );

    #include "addDictOption.H"

    argList::addNote
    (
        "test intersect of two surfaces. Writes obj file"
    );
    argList::noParallel();
    argList::noFunctionObjects();
    argList::validArgs.append("surface file");
    argList::validArgs.append("surface file");

    #include "setRootCase.H"
    #include "createTime.H"

    const word outputFile(args.executable() + ".obj");

    const fileName surf1Name(args[1]);
    triSurface surf1 = loadSurface(runTime, surf1Name)();
    Info<< surf1Name << " statistics:" << endl;
    surf1.writeStats(Info);
    Info<< endl;

    const fileName surf2Name(args[2]);
    triSurface surf2 = loadSurface(runTime, surf2Name)();
    Info<< surf2Name << " statistics:" << endl;
    surf2.writeStats(Info);
    Info<< endl;


    if (args.optionFound("debug2"))
    {
        surfaceIntersection::debug |= 2;
    }
    if (args.optionFound("debug4"))
    {
        surfaceIntersection::debug |= 4;
    }
    const bool optPrint = args.optionFound("print");

    dictionary intersectOptions;
    if (args.optionFound("dict"))
    {
        intersectOptions = IOdictionary
        (
            IOobject
            (
                args["dict"],
                runTime,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        );
    }

    intersectOptions.writeEntry("intersectOptions", Info);
    Info<< endl;

    triSurfaceSearch query1(surf1);
    triSurfaceSearch query2(surf2);
    surfaceIntersection cuts(query1, query2, intersectOptions);

    Info<<"intersection "
        << cuts.cutPoints().size() << " points "
        << cuts.cutEdges().size()  << " edges" << nl;

    if (optPrint)
    {
        Info<< "surf1-cuts: " << cuts.surf1EdgeCuts() << nl
            << "surf2-cuts: " << cuts.surf2EdgeCuts() << nl
            << "face-pairs: " << cuts.facePairToEdge() << nl
            << "edges: " << cuts.cutEdges() << nl;
    }

    word mergeOp;
    if (args.optionFound("mergePoints"))
    {
        cuts.mergePoints(args.optionRead<scalar>("mergePoints"));
        mergeOp = "mergePoints";
    }
    else if (args.optionFound("mergeEdges"))
    {
        cuts.mergeEdges();
        mergeOp = "mergeEdges";
    }

    if (!mergeOp.empty())
    {
        Info<< mergeOp << ": "
            << cuts.cutPoints().size() << " points "
            << cuts.cutEdges().size()  << " edges" << nl;

        if (optPrint)
        {
            Info<< "surf1-cuts: " << cuts.surf1EdgeCuts() << nl
                << "surf2-cuts: " << cuts.surf2EdgeCuts() << nl
                << "face-pairs: " << cuts.facePairToEdge() << nl
                << "edges: " << cuts.cutEdges() << nl;
        }
    }

    const pointField& points = cuts.cutPoints();
    const edgeList&    edges = cuts.cutEdges();

    if (points.size() || edges.size())
    {
        Info<<"write to " << outputFile << nl;

        OFstream os(outputFile);

        forAll(points, pointi)
        {
            const point& pt = points[pointi];
            os << "v " << pt.x() << ' ' << pt.y() << ' ' << pt.z() << nl;
        }

        forAll(edges, edgei)
        {
            const edge& e = edges[edgei];
            os << "l " << e.start()+1 << ' ' << e.end()+1 << nl;
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
