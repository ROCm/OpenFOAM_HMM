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

#include "STLfileFormat.H"
#include "clock.H"
#include "OSspecific.H"
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
    STLfileFormat,
    fileExtension,
    stl
);

addNamedToRunTimeSelectionTable
(
    keyedSurface,
    STLfileFormat,
    fileExtension,
    stlb
);

addNamedToMemberFunctionSelectionTable
(
    keyedSurface,
    STLfileFormat,
    write,
    fileExtension,
    stl
);

addNamedToMemberFunctionSelectionTable
(
    keyedSurface,
    STLfileFormat,
    write,
    fileExtension,
    stlb
);

addNamedToMemberFunctionSelectionTable
(
    meshedSurface,
    STLfileFormat,
    write,
    fileExtension,
    stl
);

addNamedToMemberFunctionSelectionTable
(
    meshedSurface,
    STLfileFormat,
    write,
    fileExtension,
    stlb
);

}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// check binary by getting the header and number of facets
// this seems to work better than the old token-based method
// - some programs (eg, pro-STAR) have 'solid' as the first word in
//   the binary header.
// - using wordToken can cause an abort if non-word (binary) content
//   is detected ... this is not exactly what we want.
bool Foam::fileFormats::STLfileFormat::detectBINARY
(
    const fileName& fName
)
{
    off_t fileSize = Foam::size(fName);

    IFstream ifs(fName, IOstream::BINARY);
    istream& is = ifs.stdStream();

    // Read the STL header
    char header[headerSize];
    is.read(header, headerSize);

    // Check that stream is OK, if not this may be an ASCII file
    if (!is.good())
    {
        return false;
    }

    // Read the number of triangles in the STl file
    // (note: read as int so we can check whether >2^31)
    int nTris;
    is.read(reinterpret_cast<char*>(&nTris), sizeof(unsigned int));

    // Check that stream is OK and number of triangles is positive,
    // if not this maybe an ASCII file
    //
    // Also compare the file size with that expected from the number of tris
    // If the comparison is not sensible then it may be an ASCII file
    if
    (
        !is
     || nTris < 0
     || nTris < (fileSize - headerSize)/50
     || nTris > (fileSize - headerSize)/25
    )
    {
        return false;
    }

    // looks like it might be BINARY
    return true;
}


#undef DEBUG_STLBINARY

bool Foam::fileFormats::STLfileFormat::readBINARY
(
    IFstream& ifs,
    const off_t fileSize
)
{
    istream& is = ifs.stdStream();

    // Read the STL header
    char header[headerSize];
    is.read(header, headerSize);

    // Check that stream is OK, if not this may be an ASCII file
    if (!is.good())
    {
        FatalErrorIn("fileFormats::STLfileFormat::readBINARY(Istream&)")
            << "problem reading header, perhaps file is not binary "
            << exit(FatalError);
    }

    // Read the number of triangles in the STl file
    // (note: read as int so we can check whether >2^31)
    int nTris;
    is.read(reinterpret_cast<char*>(&nTris), sizeof(unsigned int));

    // Check that stream is OK and number of triangles is positive,
    // if not this maybe an ASCII file
    //
    // Also compare the file size with that expected from the number of tris
    // If the comparison is not sensible then it may be an ASCII file
    if
    (
        !is
     || nTris < 0
     || nTris < (fileSize - headerSize)/50
     || nTris > (fileSize - headerSize)/25
    )
    {
        FatalErrorIn("fileFormats::STLfileFormat::readBINARY(Istream&)")
            << "problem reading number of triangles, perhaps file is not binary"
            << exit(FatalError);
    }

#ifdef DEBUG_STLBINARY
    Info<< "# " << nTris << " facets" << endl;
    label prevRegion = -1;
#endif


    // write directly into the lists:
    pointField& pointLst = points();
    List<face>& faceLst = faces();
    List<label>& regionLst = regions();

    pointLst.setSize(3*nTris);
    faceLst.setSize(nTris);
    regionLst.setSize(nTris);

    label maxRegionId = 0;

    label pointI = 0;
    forAll(faceLst, faceI)
    {
        face& fTri = faceLst[faceI];
        fTri.setSize(3);


        // Read an STL triangle
        STLtriangle stlTri(is);

        // Set the rawPoints to the vertices of the STL triangle
        // and set the point labels of the face

        pointLst[pointI] = stlTri.a();
        fTri[0] = pointI++;

        pointLst[pointI] = stlTri.b();
        fTri[1] = pointI++;

        pointLst[pointI] = stlTri.c();
        fTri[2] = pointI++;

        if (maxRegionId < stlTri.region())
        {
            maxRegionId = stlTri.region();
        }

        // interprete colour as a region
        regionLst[faceI] = stlTri.region();

#ifdef DEBUG_STLBINARY
        if
        (
            prevRegion != stlTri.region()
        )
        {
            if (prevRegion != -1)
            {
                Info<< "endsolid region" << prevRegion << nl;
            }
            prevRegion = stlTri.region();

            Info<< "solid region" << prevRegion << nl;
        }

        Info<< "  facet normal " << stlTri.normal() << nl;
        Info<< "    outer loop" << nl;
        Info<< "      vertex " << stlTri.a() << nl;
        Info<< "      vertex " << stlTri.b() << nl;
        Info<< "      vertex " << stlTri.c() << nl;
        Info<< "    outer loop" << nl;
        Info<< "  endfacet" << endl;
#endif
    }
#ifdef DEBUG_STLBINARY
    Info<< "endsolid region" << prevRegion << nl;
#endif

    setPatches(maxRegionId);
    stitchFaces(SMALL);

    return true;
}


inline void Foam::fileFormats::STLfileFormat::writeShell
(
    Ostream& os,
    const pointField& pointLst,
    const face& f,
    const vector& norm
)
{
    // simple triangulation about f[0].
    // better triangulation should have been done before
    const point& p0 = pointLst[f[0]];
    for (label fp1 = 1; fp1 < f.size() - 1; ++fp1)
    {
        label fp2 = (fp1 + 1) % f.size();

        const point& p1 = pointLst[f[fp1]];
        const point& p2 = pointLst[f[fp2]];

        // write STL triangle
        os  << "  facet normal "
            << norm.x() << ' ' << norm.y() << ' ' << norm.z() << nl
            << "    outer loop\n"
            << "       vertex "
            << p0.x() << ' ' << p0.y() << ' ' << p0.z() << nl
            << "       vertex "
            << p1.x() << ' ' << p1.y() << ' ' << p1.z() << nl
            << "       vertex "
            << p2.x() << ' ' << p2.y() << ' ' << p2.z() << nl
            << "    endloop\n"
            << "  endfacet" << endl;
    }
}


inline void Foam::fileFormats::STLfileFormat::writeShell
(
    ostream& os,
    const pointField& pointLst,
    const face& f,
    const vector& norm,
    const label patchI
)
{
    // simple triangulation about f[0].
    // better triangulation should have been done before
    const point& p0 = pointLst[f[0]];
    for (label fp1 = 1; fp1 < f.size() - 1; ++fp1)
    {
        label fp2 = (fp1 + 1) % f.size();

        STLtriangle stlTri
        (
            norm,
            p0,
            pointLst[f[fp1]],
            pointLst[f[fp2]],
            patchI
        );

        stlTri.write(os);
    }
}




// write sorted:
void Foam::fileFormats::STLfileFormat::writeASCII
(
    Ostream& os,
    const keyedSurface& surf
)
{
    const pointField& pointLst = surf.points();
    const List<face>& faceLst  = surf.faces();
    const vectorField& normLst = surf.faceNormals();

    labelList faceMap;
    List<surfacePatch> patchLst = surf.sortedRegions(faceMap);

    label faceIndex = 0;
    forAll(patchLst, patchI)
    {
        // Print all faces belonging to this region
        const surfacePatch& patch = patchLst[patchI];

        os << "solid " << patch.name() << endl;

        forAll(patch, patchFaceI)
        {
            const label faceI = faceMap[faceIndex++];
            writeShell(os, pointLst, faceLst[faceI], normLst[faceI]);
        }

        os << "endsolid " << patch.name() << endl;
    }
}



// write sorted:
void Foam::fileFormats::STLfileFormat::writeASCII
(
    Ostream& os,
    const meshedSurface& surf
)
{
    const pointField& pointLst = surf.points();
    const List<face>& faceLst = surf.faces();
    const List<surfacePatch>& patchLst = surf.patches();
    const vectorField& normLst = surf.faceNormals();

    // force triangulation, but just do the cheapest form possible
    label faceIndex = 0;
    forAll(patchLst, patchI)
    {
        // Print all faces belonging to this region
        const surfacePatch& patch = patchLst[patchI];

        os << "solid " << patch.name() << endl;

        forAll(patch, patchFaceI)
        {
            const label faceI = faceIndex++;
            writeShell(os, pointLst, faceLst[faceI], normLst[faceI]);
        }

        os << "endsolid " << patch.name() << endl;
    }
}


// write unsorted:
void Foam::fileFormats::STLfileFormat::writeBINARY
(
    ostream& os,
    const keyedSurface& surf
)
{
    const pointField& pointLst = surf.points();
    const List<face>&  faceLst = surf.faces();
    const List<label>& regionLst = surf.regions();
    const vectorField& normLst = surf.faceNormals();

    // Write the STL header
    string header("STL binary file", headerSize);

    // clear possible trailing junk
    for (label i = header.size(); i < headerSize; ++i)
    {
        header[i] = 0;
    }
    os.write(header.c_str(), headerSize);

    // force triangulation, but just do the cheapest form possible
    unsigned int nTris = 0;
    forAll(faceLst, faceI)
    {
        nTris += faceLst[faceI].size() - 2;
    }

    os.write(reinterpret_cast<char*>(&nTris), sizeof(unsigned int));

    // always write unsorted
    forAll(faceLst, faceI)
    {
        writeShell
        (
            os,
            pointLst,
            faceLst[faceI],
            normLst[faceI],
            regionLst[faceI]
        );
    }
}


void Foam::fileFormats::STLfileFormat::writeBINARY
(
    ostream& os,
    const meshedSurface& surf
)
{
    const pointField& pointLst = surf.points();
    const List<face>& faceLst  = surf.faces();
    const vectorField& normLst = surf.faceNormals();
    const surfacePatchList& patchLst = surf.patches();

    // Write the STL header
    string header("STL binary file", headerSize);
    os.write(header.c_str(), headerSize);

    // force triangulation, but just do the cheapest form possible

    unsigned int nTris = 0;
    forAll(faceLst, faceI)
    {
        nTris += faceLst[faceI].size() - 2;
    }

    os.write(reinterpret_cast<char*>(&nTris), sizeof(unsigned int));

    label faceIndex = 0;
    forAll(patchLst, patchI)
    {
        forAll(patchLst[patchI], patchFaceI)
        {
            writeShell
            (
                os,
                pointLst,
                faceLst[faceIndex],
                normLst[faceIndex],
                patchI
            );

            ++faceIndex;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileFormats::STLfileFormat::STLfileFormat()
:
    Foam::keyedSurface()
{}


Foam::fileFormats::STLfileFormat::STLfileFormat
(
    const fileName& fName,
    const bool
)
:
    Foam::keyedSurface()
{
    off_t fileSize = Foam::size(fName);

    // auto-detect ascii/binary
    if (detectBINARY(fName))
    {
        readBINARY(IFstream(fName, IOstream::BINARY)(), fileSize);
    }
    else
    {
        readASCII(IFstream(fName)(), fileSize);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::fileFormats::STLfileFormat::write
(
    Ostream& os,
    const keyedSurface& surf
)
{
    writeASCII(os, surf);
}


void Foam::fileFormats::STLfileFormat::write
(
    Ostream& os,
    const meshedSurface& surf
)
{
    writeASCII(os, surf);
}


void Foam::fileFormats::STLfileFormat::write
(
    const fileName& fName,
    const keyedSurface& surf
)
{
    const word ext = fName.ext();

    // handle 'stlb' as binary directly
    if (ext == "stlb")
    {
        std::ofstream ofs(fName.c_str(), std::ios::binary);
        writeBINARY(ofs, surf);
    }
    else
    {
        writeASCII(OFstream(fName)(), surf);
    }
}


void Foam::fileFormats::STLfileFormat::write
(
    const fileName& fName,
    const meshedSurface& surf
)
{
    const word ext = fName.ext();

    // handle 'stlb' as binary directly
    if (ext == "stlb")
    {
        std::ofstream ofs(fName.c_str(), std::ios::binary);
        writeBINARY(ofs, surf);
    }
    else
    {
        writeASCII(OFstream(fName)(), surf);
    }
}

// ************************************************************************* //
