/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "foamVtuSizing.H"
#include "foamVtkCore.H"
#include "polyMesh.H"
#include "cellShape.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class LabelType, class LabelType2>
void Foam::vtk::vtuSizing::populateArrays
(
    const polyMesh& mesh,
    const vtk::vtuSizing& sizing,
    UList<uint8_t>& cellTypes,
    UList<LabelType>& vertLabels,
    UList<LabelType>& vertOffset,
    UList<LabelType>& faceLabels,
    UList<LabelType>& faceOffset,
    const enum contentType output,
    UList<LabelType2>& cellMap,
    UList<LabelType2>& addPointsIds
)
{
    // STAGE 1: Verify storage sizes

    if (cellTypes.size() != sizing.nFieldCells())
    {
        FatalErrorInFunction
            << " cellTypes size=" << cellTypes.size()
            << " expected " << sizing.nFieldCells()
            << exit(FatalError);
    }

    if (cellMap.size() != sizing.nFieldCells())
    {
        FatalErrorInFunction
            << " cellMap size=" << cellMap.size()
            << " expected " << sizing.nFieldCells()
            << exit(FatalError);
    }

    if (addPointsIds.size() != sizing.nAddPoints())
    {
        FatalErrorInFunction
            << " addPointsIds size=" << addPointsIds.size()
            << " expected " << sizing.nAddPoints()
            << exit(FatalError);
    }

    // Prefix vertLabels with the size too?
    // Also use as the size of the prefixed information
    const int prefix = (output != contentType::XML) ? 1 : 0;

    switch (output)
    {
        case contentType::LEGACY:
        {
            if (vertLabels.size() != sizing.sizeLegacy())
            {
                FatalErrorInFunction
                    << " legacy size=" << vertLabels.size()
                    << " expected " << sizing.sizeLegacy()
                    << exit(FatalError);
            }
            break;
        }
        case contentType::XML:
        {
            // XML uses connectivity/offset pair.
            if
            (
                vertLabels.size()
             != sizing.sizeXml(slotType::CELLS)
            )
            {
                FatalErrorInFunction
                    << " connectivity size=" << vertLabels.size()
                    << " expected "
                    << sizing.sizeXml(slotType::CELLS)
                    << exit(FatalError);
            }

            if
            (
                vertOffset.size()
             != sizing.sizeXml(slotType::CELLS_OFFSETS)
            )
            {
                FatalErrorInFunction
                    << " offsets size=" << vertOffset.size()
                    << " expected "
                    << sizing.sizeXml(slotType::CELLS_OFFSETS)
                    << exit(FatalError);
            }

            if (sizing.nFaceLabels())
            {
                if
                (
                    faceLabels.size()
                 != sizing.sizeXml(slotType::FACES)
                )
                {
                    FatalErrorInFunction
                        << " faces size=" << faceLabels.size()
                        << " expected "
                        << sizing.sizeXml(slotType::FACES)
                        << exit(FatalError);
                }

                if
                (
                    faceOffset.size()
                 != sizing.sizeXml(slotType::FACES_OFFSETS)
                )
                {
                    FatalErrorInFunction
                        << " facesOffsets size=" << faceOffset.size()
                        << " expected "
                        << sizing.sizeXml(slotType::FACES_OFFSETS)
                        << exit(FatalError);
                }
            }
            break;
        }
        case contentType::INTERNAL:
        {
            // VTK-internal connectivity/offset pair.
            if
            (
                vertLabels.size()
             != sizing.sizeInternal(slotType::CELLS)
            )
            {
                FatalErrorInFunction
                    << " connectivity size=" << vertLabels.size()
                    << " expected "
                    << sizing.sizeInternal(slotType::CELLS)
                    << exit(FatalError);
            }

            if
            (
                vertOffset.size()
             != sizing.sizeInternal(slotType::CELLS_OFFSETS)
            )
            {
                FatalErrorInFunction
                    << " offsets size=" << vertOffset.size()
                    << " expected "
                    << sizing.sizeInternal(slotType::CELLS_OFFSETS)
                    << exit(FatalError);
            }

            if (sizing.nFaceLabels())
            {
                if
                (
                    faceLabels.size()
                 != sizing.sizeInternal(slotType::FACES)
                )
                {
                    FatalErrorInFunction
                        << " faces size=" << faceLabels.size()
                        << " expected "
                        << sizing.sizeInternal(slotType::FACES)
                        << exit(FatalError);
                }

                if
                (
                    faceOffset.size()
                 != sizing.sizeInternal(slotType::FACES_OFFSETS)
                )
                {
                    FatalErrorInFunction
                        << " facesOffsets size=" << faceOffset.size()
                        << " expected "
                        << sizing.sizeInternal(slotType::FACES_OFFSETS)
                        << exit(FatalError);
                }
            }
            break;
        }
    }


    faceOffset = -1;

    const cellModel& tet      = cellModel::ref(cellModel::TET);
    const cellModel& pyr      = cellModel::ref(cellModel::PYR);
    const cellModel& prism    = cellModel::ref(cellModel::PRISM);
    const cellModel& wedge    = cellModel::ref(cellModel::WEDGE);
    const cellModel& tetWedge = cellModel::ref(cellModel::TETWEDGE);
    const cellModel& hex      = cellModel::ref(cellModel::HEX);

    const cellShapeList& shapes = mesh.cellShapes();

    // face owner is needed to determine the face orientation
    const labelList& owner = mesh.faceOwner();

    // Unique vertex labels per polyhedral
    HashSet<label> hashUniqId(2*256);

    // Index into vertLabels, faceLabels for normal cells
    label nVertLabels = 0;
    label nFaceLabels = 0;

    // Index into vertLabels for decomposed polys
    label nVertDecomp = sizing.nVertLabels() + prefix*sizing.nCells();

    // Placement of decomposed cells
    label nCellDecomp = mesh.nCells();

    // Placement of additional point labels
    label nPointDecomp = 0;

    // Non-decomposed polyhedral are represented as a face-stream.
    // For legacy format, this stream replaces the normal connectivity
    // information. Use references to alias where the face output should land.

    UList<LabelType>& faceOutput =
    (
        output == contentType::LEGACY
      ? vertLabels
      : faceLabels
    );

    label& faceIndexer =
    (
        output == contentType::LEGACY
      ? nVertLabels
      : nFaceLabels
    );

    // ===========================================
    // STAGE 2: Rewrite in VTK form
    // During this stage, the vertOffset contains the *size* associated with
    // the per-cell vertLabels entries, and the faceOffset contains the *size*
    // associated with the per-cell faceLabels.

    forAll(shapes, celli)
    {
        const cellShape& shape = shapes[celli];
        const cellModel& model = shape.model();

        cellMap[celli] = celli;

        if (model == tet)
        {
            cellTypes[celli] = vtk::cellType::VTK_TETRA;
            if (vertOffset.size())
            {
                vertOffset[celli] = shape.size();
            }
            if (prefix)
            {
                vertLabels[nVertLabels++] = shape.size();
            }

            forAll(shape, i)
            {
                vertLabels[nVertLabels++] = shape[i];
            }
        }
        else if (model == pyr)
        {
            cellTypes[celli] = vtk::cellType::VTK_PYRAMID;
            if (vertOffset.size())
            {
                vertOffset[celli] = shape.size();
            }
            if (prefix)
            {
                vertLabels[nVertLabels++] = shape.size();
            }

            forAll(shape, i)
            {
                vertLabels[nVertLabels++] = shape[i];
            }
        }
        else if (model == hex)
        {
            cellTypes[celli] = vtk::cellType::VTK_HEXAHEDRON;
            if (vertOffset.size())
            {
                vertOffset[celli] = shape.size();
            }
            if (prefix)
            {
                vertLabels[nVertLabels++] = shape.size();
            }

            forAll(shape, i)
            {
                vertLabels[nVertLabels++] = shape[i];
            }
        }
        else if (model == prism)
        {
            cellTypes[celli] = vtk::cellType::VTK_WEDGE;
            if (vertOffset.size())
            {
                vertOffset[celli] = shape.size();
            }
            if (prefix)
            {
                vertLabels[nVertLabels++] = shape.size();
            }

            // VTK_WEDGE triangles point outwards (swap 1<->2, 4<->5)
            vertLabels[nVertLabels++] = shape[0];
            vertLabels[nVertLabels++] = shape[2];
            vertLabels[nVertLabels++] = shape[1];
            vertLabels[nVertLabels++] = shape[3];
            vertLabels[nVertLabels++] = shape[5];
            vertLabels[nVertLabels++] = shape[4];
        }
        else if (model == tetWedge && sizing.decompose())
        {
            // Treat as squeezed prism
            cellTypes[celli] = vtk::cellType::VTK_WEDGE;
            if (vertOffset.size())
            {
                vertOffset[celli] = 6;
            }
            if (prefix)
            {
                vertLabels[nVertLabels++] = 6;
            }

            vertLabels[nVertLabels++] = shape[0];
            vertLabels[nVertLabels++] = shape[2];
            vertLabels[nVertLabels++] = shape[1];
            vertLabels[nVertLabels++] = shape[3];
            vertLabels[nVertLabels++] = shape[4];
            vertLabels[nVertLabels++] = shape[3];
        }
        else if (model == wedge && sizing.decompose())
        {
            // Treat as squeezed hex
            cellTypes[celli] = vtk::cellType::VTK_HEXAHEDRON;
            if (vertOffset.size())
            {
                vertOffset[celli] = 8;
            }
            if (prefix)
            {
                vertLabels[nVertLabels++] = 8;
            }

            vertLabels[nVertLabels++] = shape[0];
            vertLabels[nVertLabels++] = shape[1];
            vertLabels[nVertLabels++] = shape[2];
            vertLabels[nVertLabels++] = shape[2];
            vertLabels[nVertLabels++] = shape[3];
            vertLabels[nVertLabels++] = shape[4];
            vertLabels[nVertLabels++] = shape[5];
            vertLabels[nVertLabels++] = shape[6];
        }
        else if (sizing.decompose())
        {
            // Polyhedral cell - decompose into tet/pyr.

            // Ensure we have the correct orientation for the base of the
            // primitive cell shape.
            // If the cell is face owner, the orientation needs to be flipped
            // to avoid defining negative cells.
            // VTK may not care, but we'll do it anyhow for safety.

            // Mapping from additional point to cell, and the new vertex from
            // the cell-centre
            const label newVertexLabel = mesh.nPoints() + nPointDecomp;

            addPointsIds[nPointDecomp++] = celli;

            // Whether to insert cell in place of original or not.
            bool first = true;

            const labelList& cFaces = mesh.cells()[celli];
            forAll(cFaces, cFaceI)
            {
                const face& f = mesh.faces()[cFaces[cFaceI]];
                const bool isOwner = (owner[cFaces[cFaceI]] == celli);

                // Count triangles/quads in decomposition
                label nTria = 0, nQuad = 0;
                f.nTrianglesQuads(mesh.points(), nTria, nQuad);

                // Do actual decomposition
                faceList faces3(nTria);
                faceList faces4(nQuad);
                nTria = 0, nQuad = 0;
                f.trianglesQuads(mesh.points(), nTria, nQuad, faces3, faces4);

                forAll(faces4, fci)
                {
                    // Quad becomes a pyramid
                    const face& quad = faces4[fci];
                    const label nShapePoints = 5;  // pyr (5 vertices)

                    label celLoc, vrtLoc;
                    if (first)
                    {
                        first = false;
                        celLoc = celli;
                        vrtLoc = nVertLabels;
                        nVertLabels += prefix + nShapePoints;
                    }
                    else
                    {
                        celLoc = nCellDecomp++;
                        vrtLoc = nVertDecomp;
                        nVertDecomp += prefix + nShapePoints;
                    }
                    cellMap[celLoc] = celli;

                    cellTypes[celLoc] = vtk::cellType::VTK_PYRAMID;
                    if (vertOffset.size())
                    {
                        vertOffset[celLoc] = nShapePoints;
                    }
                    if (prefix)
                    {
                        vertLabels[vrtLoc++] = nShapePoints;
                    }

                    // See note above about the orientation.
                    if (isOwner)
                    {
                        vertLabels[vrtLoc++] = quad[3];
                        vertLabels[vrtLoc++] = quad[2];
                        vertLabels[vrtLoc++] = quad[1];
                        vertLabels[vrtLoc++] = quad[0];
                    }
                    else
                    {
                        vertLabels[vrtLoc++] = quad[0];
                        vertLabels[vrtLoc++] = quad[1];
                        vertLabels[vrtLoc++] = quad[2];
                        vertLabels[vrtLoc++] = quad[3];
                    }

                    vertLabels[vrtLoc++] = newVertexLabel; // apex
                }

                forAll(faces3, fci)
                {
                    // Triangle becomes a tetrahedral
                    const face& tria = faces3[fci];
                    const label nShapePoints = 4;  // tet (4 vertices)

                    label celLoc, vrtLoc;
                    if (first)
                    {
                        first = false;
                        celLoc = celli;
                        vrtLoc = nVertLabels;
                        nVertLabels += prefix + nShapePoints;
                    }
                    else
                    {
                        celLoc = nCellDecomp++;
                        vrtLoc = nVertDecomp;
                        nVertDecomp += prefix + nShapePoints;
                    }
                    cellMap[celLoc] = celli;

                    cellTypes[celLoc] = vtk::cellType::VTK_TETRA;
                    if (vertOffset.size())
                    {
                        vertOffset[celLoc] = nShapePoints;
                    }
                    if (prefix)
                    {
                        vertLabels[vrtLoc++] = nShapePoints;
                    }

                    cellTypes[celLoc] = vtk::cellType::VTK_TETRA;

                    // See note above about the orientation.
                    if (isOwner)
                    {
                        vertLabels[vrtLoc++] = tria[2];
                        vertLabels[vrtLoc++] = tria[1];
                        vertLabels[vrtLoc++] = tria[0];
                    }
                    else
                    {
                        vertLabels[vrtLoc++] = tria[0];
                        vertLabels[vrtLoc++] = tria[1];
                        vertLabels[vrtLoc++] = tria[2];
                    }
                    vertLabels[vrtLoc++] = newVertexLabel; // apex
                }
            }
        }
        else
        {
            // Polyhedral cell - not decomposed
            hashUniqId.clear();  // unique node ids used (XML, INTERNAL)

            // face-stream
            //   [nFaces, nFace0Pts, id1, id2, ..., nFace1Pts, id1, id2, ...]
            cellTypes[celli] = vtk::cellType::VTK_POLYHEDRON;
            const labelList& cFaces = mesh.cells()[celli];

            const label startLabel = faceIndexer;

            if (output == contentType::LEGACY)
            {
                faceOutput[startLabel] = 0; // placeholder for size
                ++faceIndexer;
            }

            faceOutput[faceIndexer++] = cFaces.size();

            forAll(cFaces, cFaceI)
            {
                const face& f = mesh.faces()[cFaces[cFaceI]];
                const bool isOwner = (owner[cFaces[cFaceI]] == celli);

                forAll(f, fp)
                {
                    hashUniqId.insert(f[fp]);
                }

                // number of labels for this face
                faceOutput[faceIndexer++] = f.size();

                if (isOwner)
                {
                    forAll(f, fp)
                    {
                        faceOutput[faceIndexer++] = f[fp];
                    }
                }
                else
                {
                    // fairly immaterial if we reverse the list
                    // or use face::reverseFace()
                    forAllReverse(f, fp)
                    {
                        faceOutput[faceIndexer++] = f[fp];
                    }
                }
            }

            if (output == contentType::LEGACY)
            {
                // Update size for legacy face stream
                faceOutput[startLabel] = (faceIndexer - startLabel);
            }
            else
            {
                // Size for face stream
                faceOffset[celli] = (faceIndexer - startLabel);

                vertOffset[celli] = hashUniqId.size();
                if (prefix)
                {
                    vertLabels[nVertLabels++] = hashUniqId.size();
                }

                const labelList uniq = hashUniqId.sortedToc();
                forAll(uniq, i)
                {
                    vertLabels[nVertLabels++] = uniq[i];
                }
            }
        }
    }

    // ===========================================
    // STAGE 3: Adjust vertOffset for all cells
    // A second pass is needed for several reasons.
    // - Additional (decomposed) cells are placed out of sequence
    // - Internal format has the size prefixed, XML format does not.
    // - XML format expects end-offsets, Internal expects begin-offsets

    switch (output)
    {
        case contentType::LEGACY: // nothing to do
            break;

        case contentType::XML:
        {
            // No prefix, determine end offsets
            // vertOffset[0] already contains its size
            for (label i = 1; i < vertOffset.size(); ++i)
            {
                vertOffset[i] += vertOffset[i-1];
            }

            if (sizing.nFaceLabels())
            {
                // End face offsets, leaving -1 untouched
                label prev = 0;
                forAll(faceOffset, i)
                {
                    const label sz = faceOffset[i];
                    if (sz > 0)
                    {
                        prev += sz;
                        faceOffset[i] = prev;
                    }
                }
            }
            break;
        }
        case contentType::INTERNAL:
        {
            // Has prefix, determine begin offsets
            label beg = 0;
            forAll(vertOffset, i)
            {
                const label sz = vertOffset[i];
                vertOffset[i] = beg;
                beg += 1 + sz;
            }

            // Begin face offsets, leaving -1 untouched
            if (sizing.nFaceLabels())
            {
                beg = 0;
                forAll(faceOffset, i)
                {
                    const label sz = faceOffset[i];
                    if (sz > 0)
                    {
                        faceOffset[i] = beg;
                        beg += sz;
                    }
                }
            }
            break;
        }
    }
}


// ************************************************************************* //
