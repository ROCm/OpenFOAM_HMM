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
    const UnsortedMeshedSurface<Face>& surf
)
{
    const pointField& pointLst = surf.points();
    const List<Face>& faceLst  = surf.faces();
    const vectorField& normLst = surf.faceNormals();

    labelList faceMap;
    List<surfGroup> patchLst = surf.sortedRegions(faceMap);

    label faceIndex = 0;
    forAll(patchLst, patchI)
    {
        // Print all faces belonging to this region
        const surfGroup& patch = patchLst[patchI];

        os  << "solid " << patch.name() << endl;
        forAll(patch, patchFaceI)
        {
            const label faceI = faceMap[faceIndex++];
            writeShell(os, pointLst, faceLst[faceI], normLst[faceI]);
        }

        os  << "endsolid " << patch.name() << endl;
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

    // force triangulation, but just do the cheapest form possible
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

    // Write the STL header
    STLsurfaceFormatCore::writeHeaderBINARY(os);

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

    // Write the STL header
    STLsurfaceFormatCore::writeHeaderBINARY(os);

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

template<class Face>
Foam::fileFormats::STLsurfaceFormat<Face>::STLsurfaceFormat()
:
    ParentType()
{}


template<class Face>
Foam::fileFormats::STLsurfaceFormat<Face>::STLsurfaceFormat
(
    const fileName& fName
)
:
    ParentType()
{
    ThisType::read(fName);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
bool Foam::fileFormats::STLsurfaceFormat<Face>::read
(
    const fileName& fName
)
{
    ParentType::clear();

    // read in the values
    STLsurfaceFormatCore reader(fName);

    pointField&  pointLst  = ParentType::points();
    List<Face>&  faceLst   = ParentType::faces();
    List<label>& regionLst = ParentType::regions();

    // transfer
    pointLst.transfer(reader.points());
    regionLst.transfer(reader.regions());

    // assemble the faces:
    faceLst.setSize(regionLst.size());

    label ptI = 0;
    forAll(faceLst, faceI)
    {
        triFace fTri;

        fTri[0] = ptI++;
        fTri[1] = ptI++;
        fTri[2] = ptI++;

        faceLst[faceI] = fTri;
    }

    if (reader.binary())
    {
        ParentType::setPatches(reader.maxRegionId());
    }
    else
    {
        ParentType::setPatches(reader.groupToPatch());
    }

    ParentType::stitchFaces(SMALL);
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
    const fileName& fName,
    const UnsortedMeshedSurface<Face>& surf
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


template<class Face>
void Foam::fileFormats::STLsurfaceFormat<Face>::write
(
    const fileName& fName,
    const MeshedSurface<Face>& surf
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
