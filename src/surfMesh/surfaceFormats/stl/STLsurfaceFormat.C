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

#include "STLsurfaceFormat.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Face>
inline void Foam::fileFormats::STLsurfaceFormat<Face>::writeShell
(
    Ostream& os,
    const pointField& pointLst,
    const Face& f,
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
        os  << " facet normal "
            << norm.x() << ' ' << norm.y() << ' ' << norm.z() << nl
            << "  outer loop\n"
            << "   vertex " << p0.x() << ' ' << p0.y() << ' ' << p0.z() << nl
            << "   vertex " << p1.x() << ' ' << p1.y() << ' ' << p1.z() << nl
            << "   vertex " << p2.x() << ' ' << p2.y() << ' ' << p2.z() << nl
            << "  endloop\n"
            << " endfacet" << endl;
    }
}


template<class Face>
inline void Foam::fileFormats::STLsurfaceFormat<Face>::writeShell
(
    ostream& os,
    const pointField& pointLst,
    const Face& f,
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
template<class Face>
void Foam::fileFormats::STLsurfaceFormat<Face>::writeASCII
(
    Ostream& os,
    const MeshedSurface<Face>& surf
)
{
    const pointField& pointLst = surf.points();
    const List<Face>& faceLst = surf.faces();
    const List<surfGroup>& patchLst = surf.patches();
    const vectorField& normLst = surf.faceNormals();

    label faceIndex = 0;
    forAll(patchLst, patchI)
    {
        // Print all faces belonging to this region
        const surfGroup& patch = patchLst[patchI];

        os << "solid " << patch.name() << endl;
        forAll(patch, patchFaceI)
        {
            const label faceI = faceIndex++;
            writeShell(os, pointLst, faceLst[faceI], normLst[faceI]);
        }
        os << "endsolid " << patch.name() << endl;
    }
}


// write sorted:
template<class Face>
void Foam::fileFormats::STLsurfaceFormat<Face>::writeASCII
(
    Ostream& os,
    const UnsortedMeshedSurface<Face>& surf
)
{
    const pointField& pointLst = surf.points();
    const List<Face>& faceLst  = surf.faces();
    const vectorField& normLst = surf.faceNormals();

    if (surf.patches().size() == 1)
    {
        // a single region - we can skip sorting
        os << "solid " << surf.patches()[0].name() << endl;
        forAll(faceLst, faceI)
        {
            writeShell(os, pointLst, faceLst[faceI], normLst[faceI]);
        }
        os << "endsolid " << surf.patches()[0].name() << endl;
    }
   else
   {
        labelList faceMap;
        List<surfGroup> patchLst = surf.sortedRegions(faceMap);

        label faceIndex = 0;
        forAll(patchLst, patchI)
        {
            // Print all faces belonging to this region
            const surfGroup& patch = patchLst[patchI];

            os << "solid " << patch.name() << endl;
            forAll(patch, patchFaceI)
            {
                const label faceI = faceMap[faceIndex++];
                writeShell(os, pointLst, faceLst[faceI], normLst[faceI]);
            }
            os << "endsolid " << patch.name() << endl;
        }
   }
}



template<class Face>
void Foam::fileFormats::STLsurfaceFormat<Face>::writeBINARY
(
    ostream& os,
    const MeshedSurface<Face>& surf
)
{
    const pointField& pointLst = surf.points();
    const List<Face>& faceLst  = surf.faces();
    const vectorField& normLst = surf.faceNormals();
    const List<surfGroup>& patchLst = surf.patches();

    unsigned int nTris = 0;
    if (surf.isTri())
    {
        nTris = faceLst.size();
    }
    else
    {
        // count triangles for on-the-fly triangulation
        forAll(faceLst, faceI)
        {
            nTris += faceLst[faceI].size() - 2;
        }
    }

    // Write the STL header
    STLsurfaceFormatCore::writeHeaderBINARY(os, nTris);

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


// write unsorted:
template<class Face>
void Foam::fileFormats::STLsurfaceFormat<Face>::writeBINARY
(
    ostream& os,
    const UnsortedMeshedSurface<Face>& surf
)
{
    const pointField& pointLst = surf.points();
    const List<Face>&  faceLst = surf.faces();
    const List<label>& regionLst = surf.regions();
    const vectorField& normLst = surf.faceNormals();

    unsigned int nTris = 0;
    if (surf.isTri())
    {
        nTris = faceLst.size();
    }
    else
    {
        // count triangles for on-the-fly triangulation
        forAll(faceLst, faceI)
        {
            nTris += faceLst[faceI].size() - 2;
        }
    }

    // Write the STL header
    STLsurfaceFormatCore::writeHeaderBINARY(os, nTris);

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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::fileFormats::STLsurfaceFormat<Face>::STLsurfaceFormat
(
    const fileName& filename
)
{
    read(filename);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
bool Foam::fileFormats::STLsurfaceFormat<Face>::read
(
    const fileName& filename
)
{
    this->clear();

    // read in the values
    STLsurfaceFormatCore reader(filename);

    // transfer points
    this->storedPoints().transfer(reader.points());

    // retrieve the original region information
    List<word>  names(xferMove(reader.names()));
    List<label> sizes(xferMove(reader.sizes()));
    List<label> regions(xferMove(reader.regions()));

    // generate the (sorted) faces
    List<Face> faceLst(regions.size());

    if (reader.sorted())
    {
        // already sorted - generate directly
        forAll(faceLst, faceI)
        {
            const label startPt = 3*faceI;
            faceLst[faceI] = triFace(startPt, startPt+1, startPt+2);
        }
    }
    else
    {
        // unsorted - determine the sorted order:
        // avoid SortableList since we discard the main list anyhow
        List<label> faceMap;
        sortedOrder(regions, faceMap);

        // generate sorted faces
        forAll(faceMap, faceI)
        {
            const label startPt = 3*faceMap[faceI];
            faceLst[faceI] = triFace(startPt, startPt+1, startPt+2);
        }
    }
    regions.clear();

    // transfer:
    this->storedFaces().transfer(faceLst);

    if (names.size())
    {
        this->addPatches(sizes, names);
    }
    else
    {
        this->addPatches(sizes);
    }

    this->stitchFaces(SMALL);
    return true;
}


template<class Face>
void Foam::fileFormats::STLsurfaceFormat<Face>::write
(
    Ostream& os,
    const UnsortedMeshedSurface<Face>& surf
)
{
    writeASCII(os, surf);
}


template<class Face>
void Foam::fileFormats::STLsurfaceFormat<Face>::write
(
    Ostream& os,
    const MeshedSurface<Face>& surf
)
{
    writeASCII(os, surf);
}


template<class Face>
void Foam::fileFormats::STLsurfaceFormat<Face>::write
(
    const fileName& filename,
    const UnsortedMeshedSurface<Face>& surf
)
{
    word ext = filename.ext();

    // handle 'stlb' as binary directly
    if (ext == "stlb")
    {
        std::ofstream ofs(filename.c_str(), std::ios::binary);
        writeBINARY(ofs, surf);
    }
    else
    {
        writeASCII(OFstream(filename)(), surf);
    }
}


template<class Face>
void Foam::fileFormats::STLsurfaceFormat<Face>::write
(
    const fileName& filename,
    const MeshedSurface<Face>& surf
)
{
    const word ext = filename.ext();

    // handle 'stlb' as binary directly
    if (ext == "stlb")
    {
        std::ofstream ofs(filename.c_str(), std::ios::binary);
        writeBINARY(ofs, surf);
    }
    else
    {
        writeASCII(OFstream(filename)(), surf);
    }
}

// ************************************************************************* //
