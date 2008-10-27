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

#include "GTSfileFormat.H"
#include "clock.H"
#include "IFstream.H"
#include "IStringStream.H"
#include "tensor.H"
#include "primitivePatch.H"
#include "addToRunTimeSelectionTable.H"
#include "addToMemberFunctionSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //
namespace Foam
{
namespace fileFormats
{

addNamedToRunTimeSelectionTable
(
    keyedSurface,
    GTSfileFormat,
    fileExtension,
    gts
);

addNamedToMemberFunctionSelectionTable
(
    keyedSurface,
    GTSfileFormat,
    write,
    fileExtension,
    gts
);

addNamedToMemberFunctionSelectionTable
(
    meshedSurface,
    GTSfileFormat,
    write,
    fileExtension,
    gts
);

}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileFormats::GTSfileFormat::GTSfileFormat()
:
    Foam::keyedSurface()
{}


Foam::fileFormats::GTSfileFormat::GTSfileFormat
(
    const fileName& fName,
    const bool triangulate
)
:
    Foam::keyedSurface()
{
    IFstream is(fName);

    if (!is.good())
    {
        FatalErrorIn
        (
            "fileFormats::GTSfileFormat::GTSfileFormat"
            "(const fileName&)"
        )
            << "Cannot read file " << fName
            << exit(FatalError);
    }

    // Read header
    string line = getLineNoComment(is);

    label nPoints, nEdges, nElems;
    {
        IStringStream lineStream(line);
        lineStream
            >> nPoints
            >> nEdges
            >> nElems;
    }

    pointField pointLst(nPoints);

    // Read points
    forAll(pointLst, pointI)
    {
        scalar x, y, z;
        line = getLineNoComment(is);
        {
            IStringStream lineStream(line);
            lineStream
                >> x >> y >> z;
        }

        pointLst[pointI] = point(x, y, z);
    }

    // Read edges (Foam indexing)
    edgeList edges(nEdges);
    forAll(edges, edgei)
    {
        label beg, end;
        line = getLineNoComment(is);
        {
            IStringStream lineStream(line);
            lineStream
                >> beg >> end;
        }
        edges[edgei] = edge(beg - 1, end - 1);
    }


    // Read triangles. Convert references to edges into pointlabels
    List<keyedFace> faceLst(nElems);

    label maxPatch = 0;

    forAll(faceLst, faceI)
    {
        label e0Label, e1Label, e2Label;
        label region = 0;

        line = getLineNoComment(is);
        {
            IStringStream lineStream(line);
            lineStream
                >> e0Label >> e1Label >> e2Label;

            // Optional region number: read first, then check state on stream
            if (lineStream)
            {
                label num;
                lineStream >> num;
                if (!lineStream.bad())
                {
                    region = num;
                    if (maxPatch < region)
                    {
                        maxPatch = region;
                    }
                }
            }
        }

        // Determine ordering of edges e0, e1
        //  common: common vertex, shared by e0 and e1
        //  e0Far:  vertex on e0 which is not common
        //  e1Far:  vertex on e1 which is not common
        const edge& e0 = edges[e0Label - 1];
        const edge& e1 = edges[e1Label - 1];
        const edge& e2 = edges[e2Label - 1];

        label common01 = e0.commonVertex(e1);
        if (common01 == -1)
        {
            FatalErrorIn
            (
                "fileFormats::GTSfileFormat::GTSfileFormat"
                "(const fileName&)"
            )
                << "Edges 0 and 1 of triangle " << faceI
                << " do not share a point.\n"
                << "    edge0:" << e0 << nl
                << "    edge1:" << e1
                << exit(FatalError);
        }

        label e0Far = e0.otherVertex(common01);
        label e1Far = e1.otherVertex(common01);

        label common12 = e1.commonVertex(e2);
        if (common12 == -1)
        {
            FatalErrorIn
            (
                "fileFormats::GTSfileFormat::GTSfileFormat"
                "(const fileName&)"
            )
                << "Edges 1 and 2 of triangle " << faceI
                << " do not share a point.\n"
                << "    edge1:" << e1 << nl
                << "    edge2:" << e2
                << exit(FatalError);
        }
        label e2Far = e2.otherVertex(common12);

        // Does edge2 sit between edge1 and 0?
        if (common12 != e1Far || e2Far != e0Far)
        {
            FatalErrorIn
            (
                "fileFormats::GTSfileFormat::GTSfileFormat"
                "(const fileName&)"
            )
                << "Edges of triangle " << faceI
                << " reference more than three points.\n"
                << "    edge0:" << e0 << nl
                << "    edge1:" << e1 << nl
                << "    edge2:" << e2 << nl
                << exit(FatalError);
        }

        keyedFace lface(face(3), region);

        lface[0] = e0Far;
        lface[1] = common01;
        lface[2] = e1Far;

        faceLst[faceI] = lface;
    }

    // transfer to normal lists
    points().transfer(pointLst);
    faces().transfer(faceLst);
    setPatches(maxPatch);
    // stitchFaces(SMALL);
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //
// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fileFormats::GTSfileFormat::write
(
    Ostream& os,
    const keyedSurface& surf
)
{
    const pointField& pointLst = surf.points();
    const List<keyedFace>& faceLst = surf.faces();

    // It is too annoying to triangulate on-the-fly
    // just issue a warning and get out
    label nNonTris = 0;
    forAll(faceLst, faceI)
    {
        if (faceLst[faceI].size() > 3)
        {
            ++nNonTris;
        }
    }

    if (nNonTris)
    {
        FatalErrorIn
        (
            "fileFormats::GTSfileFormat::write"
            "(Ostream&, const keyedSurface&)"
        )
            << "Surface has " << nNonTris << "/" << faceLst.size()
            << " non-triangulated faces - not writing!" << endl;
        return;
    }

    labelList faceMap;
    List<surfacePatch> patchLst = surf.sortedRegions(faceMap);


    // Write header, print patch names as comment
    os  << "# GTS file" << nl
        << "# Regions:" << nl;

    forAll(patchLst, patchI)
    {
        os  << "#     " << patchI << "    "
            << patchLst[patchI].name() << nl;
    }
    os  << "#" << endl;


    os  << "# nPoints  nEdges  nTriangles" << nl
        << pointLst.size() << ' ' << surf.nEdges() << ' '
        << surf.size() << endl;


    // Write vertex coords
    forAll(pointLst, pointI)
    {
        os  << pointLst[pointI].x() << ' '
            << pointLst[pointI].y() << ' '
            << pointLst[pointI].z() << endl;
    }


    // Write edges.
    // Note: edges are in local point labels so convert
    const edgeList& es = surf.edges();
    const labelList& meshPts = surf.meshPoints();

    forAll(es, edgeI)
    {
        os  << meshPts[es[edgeI].start()] + 1 << ' '
            << meshPts[es[edgeI].end()] + 1 << endl;
    }


    // Write faces in terms of edges.
    const labelListList& faceEs = surf.faceEdges();

    label faceIndex = 0;
    forAll(patchLst, patchI)
    {
        forAll(patchLst[patchI], patchFaceI)
        {
            const labelList& fEdges = faceEs[faceMap[faceIndex++]];

            os  << fEdges[0] + 1 << ' '
                << fEdges[1] + 1 << ' '
                << fEdges[2] + 1 << ' '
                << patchI << endl;
        }
    }
}


void Foam::fileFormats::GTSfileFormat::write
(
    Ostream& os,
    const meshedSurface& surf
)
{
    const pointField& pointLst = surf.points();
    const List<surfacePatch>& patchLst = surf.patches();
    const List<face>& faceLst = surf.faces();

    // It is too annoying to triangulate on-the-fly
    // just issue a warning and get out
    label nNonTris = 0;
    forAll(faceLst, faceI)
    {
        if (faceLst[faceI].size() > 3)
        {
            ++nNonTris;
        }
    }

    if (nNonTris)
    {
        WarningIn
        (
            "fileFormats::GTSfileFormat::write"
            "(Ostream&, const meshedSurface&)"
        )
            << "Surface has " << nNonTris << "/" << faceLst.size()
            << " non-triangulated faces - not writing!" << endl;
        
        return;
    }

    // Write header, print patch names as comment
    os  << "# GTS file" << nl
        << "# Regions:" << nl;

    forAll(patchLst, patchI)
    {
        os  << "#     " << patchI << "    "
            << patchLst[patchI].name() << nl;
    }
    os  << "#" << endl;

    os  << "# nPoints  nEdges  nTriangles" << nl
        << pointLst.size() << ' ' << surf.nEdges() << ' '
        << surf.size() << endl;


    // Write vertex coords
    forAll(pointLst, pointI)
    {
        os  << pointLst[pointI].x() << ' '
            << pointLst[pointI].y() << ' '
            << pointLst[pointI].z() << endl;
    }


    // Write edges.
    // Note: edges are in local point labels so convert
    const edgeList& es = surf.edges();
    const labelList& meshPts = surf.meshPoints();

    forAll(es, edgei)
    {
        os  << meshPts[es[edgei].start()] + 1 << ' '
            << meshPts[es[edgei].end()] + 1 << endl;
    }


    // Write faces in terms of edges.
    const labelListList& faceEs = surf.faceEdges();

    label faceIndex = 0;
    forAll(patchLst, patchI)
    {
        const surfacePatch& patch = patchLst[patchI];

        forAll(patch, patchFaceI)
        {
            const labelList& fEdges = faceEs[faceIndex++];

            os  << fEdges[0] + 1 << ' '
                << fEdges[1] + 1 << ' '
                << fEdges[2] + 1 << ' '
                << patchI << endl;
        }
    }
}

// ************************************************************************* //
