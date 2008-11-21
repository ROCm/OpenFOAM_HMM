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

#include "STARCDsurfaceFormat.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Face>
inline void Foam::fileFormats::STARCDsurfaceFormat<Face>::writeShell
(
    Ostream& os,
    const Face& f,
    const label cellId,
    const label cellTableId
)
{
    os
        << cellId                    // includes 1 offset
        << " " << starcdShellShape_  // 3(shell) shape
        << " " << f.size()
        << " " << cellTableId
        << " " << starcdShellType_;  // 4(shell)

    // primitives have <= 8 vertices, but prevent overrun anyhow
    // indent following lines for ease of reading
    label count = 0;
    forAll(f, fp)
    {
        if ((count % 8) == 0)
        {
            os
                << nl
                << "  " << cellId;
        }
        os  << " " << f[fp] + 1;
        count++;
    }
    os  << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::fileFormats::STARCDsurfaceFormat<Face>::STARCDsurfaceFormat
(
    const fileName& filename
)
{
    read(filename);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
bool Foam::fileFormats::STARCDsurfaceFormat<Face>::read
(
    const fileName& filename
)
{
    const bool mustTriangulate = this->isTri();
    this->clear();

    fileName baseName = filename.lessExt();

    // STAR-CD index of points
    List<label> pointId;

    // read points from .vrt file
    readPoints
    (
        IFstream(baseName + ".vrt")(),
        this->storedPoints(),
        pointId
    );

    // Build inverse mapping (STAR-CD pointId -> index)
    Map<label> mapPointId(2*pointId.size());
    forAll(pointId, i)
    {
        mapPointId.insert(pointId[i], i);
    }
    pointId.clear();

    //
    // read .cel file
    // ~~~~~~~~~~~~~~
    IFstream is(baseName + ".cel");
    if (!is.good())
    {
        FatalErrorIn
        (
            "fileFormats::STARCDsurfaceFormat::read(const fileName&)"
        )
            << "Cannot read file " << is.name()
            << exit(FatalError);
    }

    readHeader(is, "PROSTAR_CELL");

    DynamicList<Face>  dynFaces;
    DynamicList<label> dynRegions;
    DynamicList<word>  dynNames;
    DynamicList<label> dynSizes;
    Map<label> lookup;

    // assume the cellTableIds are not intermixed
    bool sorted = true;
    label regionI = 0;

    label lineLabel, shapeId, nLabels, cellTableId, typeId;
    labelList starLabels(64);

    while ((is >> lineLabel).good())
    {
        is >> shapeId >> nLabels >> cellTableId >> typeId;

        if (nLabels > starLabels.size())
        {
            starLabels.setSize(nLabels);
        }
        starLabels = -1;

        // read indices - max 8 per line
        for (label i = 0; i < nLabels; ++i)
        {
            if ((i % 8) == 0)
            {
               is >> lineLabel;
            }
            is >> starLabels[i];
        }

        if (typeId == starcdShellType_)
        {
            // Convert groupID into patchID
            Map<label>::const_iterator fnd = lookup.find(cellTableId);
            if (fnd != lookup.end())
            {
                if (regionI != fnd())
                {
                    // cellTableIds are intermixed
                    sorted = false;
                }
                regionI = fnd();
            }
            else
            {
                regionI = dynSizes.size();
                lookup.insert(cellTableId, regionI);
                dynNames.append(word("cellTable_") + ::Foam::name(regionI));
                dynSizes.append(0);
            }

            SubList<label> vertices(starLabels, nLabels);

            // convert orig vertex id to point label
            forAll(vertices, i)
            {
                vertices[i] = mapPointId[vertices[i]];
            }

            if (mustTriangulate && nLabels > 3)
            {
                face f(vertices);

                faceList triFaces(f.nTriangles(this->points()));
                label nTri = 0;
                f.triangles(this->points(), nTri, triFaces);

                forAll(triFaces, faceI)
                {
                    // a triangular face, but not yet a triFace
                    dynFaces.append
                    (
                        triFace
                        (
                            static_cast<UList<label>&>(triFaces[faceI])
                        )
                    );
                    dynRegions.append(regionI);
                    dynSizes[regionI]++;
                }
            }
            else
            {
                dynFaces.append(Face(vertices));
                dynRegions.append(regionI);
                dynSizes[regionI]++;
            }
        }
    }
    mapPointId.clear();

    sortFacesAndStore
    (
        xferMoveTo<List<Face> >(dynFaces),
        xferMoveTo<List<label> >(dynRegions),
        sorted
    );

    // add patches, culling empty groups
    this->addPatches(dynSizes, dynNames, true);
    return true;
}


template<class Face>
void Foam::fileFormats::STARCDsurfaceFormat<Face>::write
(
    const fileName& filename,
    const MeshedSurface<Face>& surf
)
{
    fileName baseName = filename.lessExt();

    writePoints(OFstream(baseName + ".vrt")(), surf.points());
    OFstream os(baseName + ".cel");
    writeHeader(os, "CELL");

    const List<Face>& faceLst = surf.faces();
    const List<surfGroup>& patchLst = surf.patches();

    label faceIndex = 0;
    forAll(patchLst, patchI)
    {
        const surfGroup& patch = patchLst[patchI];

        forAll(patch, patchFaceI)
        {
            const Face& f = faceLst[faceIndex++];
            writeShell(os, f, faceIndex, patchI + 1);
        }
    }

    // write simple .inp file
    writeCase
    (
        OFstream(baseName + ".inp")(),
        surf.points(),
        surf.size(),
        patchLst
    );
}


template<class Face>
void Foam::fileFormats::STARCDsurfaceFormat<Face>::write
(
    const fileName& filename,
    const UnsortedMeshedSurface<Face>& surf
)
{
    fileName baseName = filename.lessExt();

    writePoints(OFstream(baseName + ".vrt")(), surf.points());

    OFstream os(baseName + ".cel");
    writeHeader(os, "CELL");

    const List<Face>& faceLst = surf.faces();
    labelList faceMap;
    List<surfGroup> patchLst = surf.sortedRegions(faceMap);

    label faceIndex = 0;
    forAll(patchLst, patchI)
    {
        const surfGroup& patch = patchLst[patchI];

        forAll(patch, patchFaceI)
        {
            const Face& f = faceLst[faceMap[faceIndex++]];
            writeShell(os, f, faceIndex, patchI + 1);
        }
    }

    // write simple .inp file
    writeCase
    (
        OFstream(baseName + ".inp")(),
        surf.points(),
        surf.size(),
        patchLst
    );
}

// ************************************************************************* //
