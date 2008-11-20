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
    const fileName& fName
)
{
    read(fName);
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
bool Foam::fileFormats::STARCDsurfaceFormat<Face>::read
(
    const fileName& fName
)
{
    const bool mustTriangulate = this->isTri();
    this->clear();

    fileName baseName = fName.lessExt();
    autoPtr<IFstream> isPtr;


    DynamicList<point> pointLst;
    // STAR-CD index of points
    DynamicList<label> pointId;

    //
    // read .vrt file
    // ~~~~~~~~~~~~~~
    isPtr.reset(new IFstream(baseName + ".vrt"));

    if (!isPtr().good())
    {
        FatalErrorIn
        (
            "fileFormats::STARCDsurfaceFormat::read(const fileName&)"
        )
            << "Cannot read file " << (baseName + ".vrt")
            << exit(FatalError);
    }

    readHeader(isPtr(), "PROSTAR_VERTEX");

    label lineLabel;

    while ((isPtr() >> lineLabel).good())
    {
        pointId.append(lineLabel);
        scalar x, y, z;

        isPtr() >> x >> y >> z;

        pointLst.append(point(x, y, z));
    }

    // transfer to normal lists
    this->storedPoints().transfer(pointLst);

    // Build inverse mapping (index to point)
    pointId.shrink();
    Map<label> mapToFoamPointId(2*pointId.size());
    forAll(pointId, i)
    {
        mapToFoamPointId.insert(pointId[i], i);
    }
    pointId.clear();


    DynamicList<Face>  faceLst;
    DynamicList<label> regionLst;

    // From face cellTableId to patchId
    Map<label> cellTableToPatchId;
    label nPatches = 0;


    //
    // read .cel file
    // ~~~~~~~~~~~~~~
    isPtr.reset(new IFstream(baseName + ".cel"));

    if (!isPtr().good())
    {
        FatalErrorIn
        (
            "fileFormats::STARCDsurfaceFormat::read(const fileName&)"
        )
            << "Cannot read file " << (baseName + ".cel")
            << exit(FatalError);
    }

    readHeader(isPtr(), "PROSTAR_CELL");

    label shapeId, nLabels, cellTableId, typeId;
    labelList starLabels(64);

    while ((isPtr() >> lineLabel).good())
    {
        isPtr() >> shapeId >> nLabels >> cellTableId >> typeId;

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
               isPtr() >> lineLabel;
            }
            isPtr() >> starLabels[i];
        }

        if (typeId == starcdShellType_)
        {
            // Convert groupID into patchID
            Map<label>::const_iterator iter =
                cellTableToPatchId.find(cellTableId);

            label patchI;
            if (iter == cellTableToPatchId.end())
            {
                patchI = nPatches++;

                cellTableToPatchId.insert(cellTableId, patchI);
            }
            else
            {
                patchI = iter();
            }


            // convert orig vertex id to point label
            for (label i=0; i < nLabels; ++i)
            {
                starLabels[i] = mapToFoamPointId[starLabels[i]];
            }

            if (mustTriangulate && nLabels > 3)
            {
                face f
                (
                    SubList<label>(starLabels, nLabels)
                );

                faceList triFaces(f.nTriangles(this->points()));
                label nTri = 0;
                f.triangles(this->points(), nTri, triFaces);

                forAll(triFaces, faceI)
                {
                    // a triangle, but not yet a triFace
                    faceLst.append
                    (
                        triFace
                        (
                            static_cast<UList<label>&>(triFaces[faceI])
                        )
                    );
                    regionLst.append(patchI);
                }
            }
            else
            {
                faceLst.append
                (
                    Face(SubList<label>(starLabels, nLabels))
                );
                regionLst.append(patchI);
            }
        }
    }

    mapToFoamPointId.clear();

    // convert cellTable_N patchId => name
    Map<word> regionNames;
    forAllConstIter(Map<label>, cellTableToPatchId, iter)
    {
        regionNames.insert
        (
            iter(),
            "cellTable_" + Foam::name(iter.key())
        );
    }

    // transfer to normal lists
    this->storedFaces().transfer(faceLst);
    this->storedRegions().transfer(regionLst);

    this->setPatches(regionNames);

    return true;
}


template<class Face>
void Foam::fileFormats::STARCDsurfaceFormat<Face>::write
(
    const fileName& fName,
    const UnsortedMeshedSurface<Face>& surf
)
{
    fileName baseName = fName.lessExt();

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


template<class Face>
void Foam::fileFormats::STARCDsurfaceFormat<Face>::write
(
    const fileName& fName,
    const MeshedSurface<Face>& surf
)
{
    fileName baseName = fName.lessExt();

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

// ************************************************************************* //
