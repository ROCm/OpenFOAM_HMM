/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "ensightMeshReader.H"
#include "cellModel.H"
#include "ensightReadFile.H"
#include "matchPoints.H"
#include "mergePoints.H"
#include "ListListOps.H"
#include "stringOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{
    defineTypeNameAndDebug(ensightMeshReader, 0);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::face& Foam::fileFormats::ensightMeshReader::rotateFace
(
    const face& f,
    face& rotatedFace
) const
{
    label fp = findMin(f);

    rotatedFace.setSize(f.size());
    forAll(rotatedFace, i)
    {
        rotatedFace[i] = f[fp];
        fp = f.fcIndex(fp);
    }
    return rotatedFace;
}


void Foam::fileFormats::ensightMeshReader::readVerts
(
    ensightReadFile& is,
    const label nVerts,
    const Map<label>& nodeIdToPoints,
    DynamicList<label>& verts
) const
{
    verts.clear();
    for (label i = 0; i < nVerts; i++)
    {
        label verti;
        is.read(verti);
        //if (nodeIdToPoints.size())
        //{
        //    verts.append(nodeIdToPoints[verti]);
        //}
        //else
        {
            verts.append(verti-1);
        }
    }
}


void Foam::fileFormats::ensightMeshReader::readIDs
(
    ensightReadFile& is,
    const bool doRead,
    const label nShapes,
    labelList& foamToElem,
    Map<label>& elemToFoam
) const
{
    const label sz = foamToElem.size();
    foamToElem.resize(sz+nShapes);
    if (doRead)
    {
        elemToFoam.resize(sz+nShapes);
        for (label shapei = 0; shapei < nShapes; shapei++)
        {
            label elemi;
            is.read(elemi);
            foamToElem[sz+shapei] = elemi;
            elemToFoam.insert(elemi, sz+shapei);
        }
    }
    else
    {
        for (label shapei = 0; shapei < nShapes; shapei++)
        {
            foamToElem[sz+shapei] = sz+shapei;
        }
    }
}


void Foam::fileFormats::ensightMeshReader::setHandedness
(
    const cellModel& model,
    DynamicList<label>& verts,
    const pointField& points
) const
{
//    // From plot3dToFoam/hexBlock.C
//    const vector x(points[verts[1]]-points[verts[0]]);
//    const scalar xMag(mag(x));
//
//    const vector y(points[verts[3]]-points[verts[0]]);
//    const scalar yMag(mag(y));
//
//    const vector z(points[verts[4]]-points[verts[0]]);
//    const scalar zMag(mag(z));
//
//    if (xMag > SMALL && yMag > SMALL && zMag > SMALL)
//    {
//        if (((x ^ y) & z) < 0)
//        {
//            // Flipped hex
//            Swap(verts[0], verts[4]);
//            Swap(verts[1], verts[5]);
//            Swap(verts[2], verts[6]);
//            Swap(verts[3], verts[7]);
//        }
//    }

    if (model.mag(verts, points) < 0)
    {
        if (verts.size() == 8)
        {
            // Flipped hex
            Swap(verts[0], verts[4]);
            Swap(verts[1], verts[5]);
            Swap(verts[2], verts[6]);
            Swap(verts[3], verts[7]);
        }
        else if (verts.size() == 4)
        {
            // Flipped tet. Change orientation of base
            Swap(verts[0], verts[1]);
        }
        else if (verts.size() == 5)
        {
            // Flipped pyr. Change orientation of base
            Swap(verts[1], verts[3]);
        }
        else if (verts.size() == 6)
        {
            // Flipped prism.
            Swap(verts[0], verts[3]);
            Swap(verts[1], verts[4]);
            Swap(verts[2], verts[5]);
        }
    }
}


bool Foam::fileFormats::ensightMeshReader::readGoldPart
(
    ensightReadFile& is,
    const bool read_node_ids,
    const bool read_elem_ids,

    pointField& points,
    labelList& pointToNodeIds,
    Map<label>& nodeIdToPoints,

    // 3D-elems : cells (cell-to-faces)
    faceListList& cells,
    labelList& cellToElemIds,
    Map<label>& elemIdToCells,

    // 2D-elems : faces
    faceList& faces,
    labelList& faceToElemIDs,
    Map<label>& elemIdToFaces
) const
{
    //- Read a single part. Return true if end-of-file reached. Return false
    //  if reaching next 'part'.


    // Work
    DynamicList<label> verts; 

    string line;
    while (is.good())
    {
        do
        {
            is.readKeyword(line);
        }
        while (line.empty() && is.good());

        const auto split = stringOps::splitSpace(line);

        if (split.size() == 0)
        {
            continue;
        }

        if (split[0] == "part")
        {
            return false;
        }
        else if (split[0] == "node_ids")
        {
            const label nPoints = points.size();
            // Ignore for now
            for (label i = 0; i < nPoints; i++)
            {
                label index;
                is.read(index);
            }
        }
        else if (split[0] == "coordinates")
        {
            label nPoints;
            is.read(nPoints);

            Pout<< indent << "coordinates " << nPoints
                << " starting at line " << is.lineNumber()
                << " position " << is.stdStream().tellg() << endl;

            readIDs
            (
                is,
                read_node_ids,
                nPoints,
                pointToNodeIds,
                nodeIdToPoints
            );

            points.setSize(nPoints);

            for (label pointi = 0; pointi < nPoints; pointi++)
            {
                is.read(points[pointi].x());
            }
            for (label pointi = 0; pointi < nPoints; pointi++)
            {
                is.read(points[pointi].y());
            }
            for (label pointi = 0; pointi < nPoints; pointi++)
            {
                is.read(points[pointi].z());
            }
        }
        else if (split[0] == "tetra4")
        {
            label nShapes;
            is.read(nShapes);

            Pout<< indent<< "tetra4 " << nShapes
                << " starting at line " << is.lineNumber()
                << " position " << is.stdStream().tellg() << endl;

            const label celli = cells.size();
            cells.resize(celli+nShapes);

            readIDs
            (
                is,
                read_elem_ids,
                nShapes,
                cellToElemIds,
                elemIdToCells
            );

            const auto& model = *cellModel::ptr(cellModel::TET);
            for (label shapei = 0; shapei < nShapes; shapei++)
            {
                readVerts(is, 4, nodeIdToPoints, verts);
                if (setHandedness_)
                {
                    setHandedness(model, verts, points);
                }
                const cellShape cellVerts(model, verts);
                cells[celli+shapei] = cellVerts.faces();
            }
        }
        else if (split[0] == "pyramid5")
        {
            label nShapes;
            is.read(nShapes);

            Pout<< indent<< "pyramid5 " << nShapes
                << " starting at line " << is.lineNumber()
                << " position " << is.stdStream().tellg() << endl;

            const label celli = cells.size();
            cells.resize(celli+nShapes);

            readIDs
            (
                is,
                read_elem_ids,
                nShapes,
                cellToElemIds,
                elemIdToCells
            );

            const auto& model = *cellModel::ptr(cellModel::PYR);
            for (label shapei = 0; shapei < nShapes; shapei++)
            {
                readVerts(is, 5, nodeIdToPoints, verts);
                if (setHandedness_)
                {
                    setHandedness(model, verts, points);
                }
                const cellShape cellVerts(model, verts);
                cells[celli+shapei] = cellVerts.faces();
            }
        }
        else if (split[0] == "penta6")
        {
            label nShapes;
            is.read(nShapes);

            Pout<< indent<< "penta6 " << nShapes
                << " starting at line " << is.lineNumber()
                << " position " << is.stdStream().tellg() << endl;

            const label celli = cells.size();
            cells.resize(celli+nShapes);

            readIDs
            (
                is,
                read_elem_ids,
                nShapes,
                cellToElemIds,
                elemIdToCells
            );

            const auto& model = *cellModel::ptr(cellModel::PRISM);
            for (label shapei = 0; shapei < nShapes; shapei++)
            {
                readVerts(is, 6, nodeIdToPoints, verts);
                if (setHandedness_)
                {
                    setHandedness(model, verts, points);
                }
                const cellShape cellVerts(model, verts);
                cells[celli+shapei] = cellVerts.faces();
            }
        }
        else if (split[0] == "hexa8")
        {
            label nShapes;
            is.read(nShapes);

            Pout<< indent<< "hexa8 " << nShapes
                << " starting at line " << is.lineNumber()
                << " position " << is.stdStream().tellg() << endl;

            const label celli = cells.size();
            cells.resize(celli+nShapes);

            readIDs
            (
                is,
                read_elem_ids,
                nShapes,
                cellToElemIds,
                elemIdToCells
            );

            const auto& model = *cellModel::ptr(cellModel::HEX);
            for (label shapei = 0; shapei < nShapes; shapei++)
            {
                readVerts(is, 8, nodeIdToPoints, verts);
                if (setHandedness_)
                {
                    setHandedness(model, verts, points);
                }
                const cellShape cellVerts(model, verts);
                cells[celli+shapei] = cellVerts.faces();
            }
        }
        else if (split[0] == "nfaced")
        {
            label nShapes;
            is.read(nShapes);

            Pout<< indent<< "nfaced " << nShapes
                << " starting at line " << is.lineNumber()
                << " position " << is.stdStream().tellg() << endl;

            const label celli = cells.size();
            cells.resize(celli+nShapes);

            readIDs
            (
                is,
                read_elem_ids,
                nShapes,
                cellToElemIds,
                elemIdToCells
            );

            for (label shapei = 0; shapei < nShapes; shapei++)
            {
                label nFaces;
                is.read(nFaces);
                faceList& cellFaces = cells[celli+shapei];
                cellFaces.setSize(nFaces);
            }

            for (label shapei = 0; shapei < nShapes; shapei++)
            {
                faceList& cellFaces = cells[celli+shapei];
                forAll(cellFaces, cellFacei)
                {
                    label nVerts;
                    is.read(nVerts);
                    cellFaces[cellFacei].setSize(nVerts);
                }
            }

            for (label shapei = 0; shapei < nShapes; shapei++)
            {
                faceList& cellFaces = cells[celli+shapei];
                forAll(cellFaces, cellFacei)
                {
                    face& f = cellFaces[cellFacei];
                    readVerts(is, f.size(), nodeIdToPoints, verts);
                    f.labelList::operator=(verts);

                    forAll(f, fp)
                    {
                        if (f[fp] < 0 || f[fp] >= points.size())
                        {
                            FatalErrorInFunction<< "Face:" << shapei
                                << " verts:" << f
                                << " indexes outside points:" << points.size()
                                << exit(FatalError);
                        }
                    }
                }
            }
        }
        else if (split[0] == "tria3")
        {
            label nShapes;
            is.read(nShapes);

            Pout<< indent << "tria3 " << nShapes
                << " starting at line " << is.lineNumber()
                << " position " << is.stdStream().tellg() << endl;

            const label facei = faces.size();

            readIDs
            (
                is,
                read_elem_ids,
                nShapes,
                faceToElemIDs,
                elemIdToFaces
            );

            faces.setSize(facei+nShapes);

            for (label shapei = 0; shapei < nShapes; shapei++)
            {
                auto& f = faces[facei+shapei];
                f.setSize(3);
                readVerts(is, f.size(), nodeIdToPoints, verts);
                f.labelList::operator=(verts);
            }
        }
        else if (split[0] == "quad4")
        {
            label nShapes;
            is.read(nShapes);

            Pout<< indent << "quad4 " << nShapes
                << " starting at line " << is.lineNumber()
                << " position " << is.stdStream().tellg() << endl;

            const label facei = faces.size();

            readIDs
            (
                is,
                read_elem_ids,
                nShapes,
                faceToElemIDs,
                elemIdToFaces
            );

            faces.setSize(facei+nShapes);

            for (label shapei = 0; shapei < nShapes; shapei++)
            {
                auto& f = faces[facei+shapei];
                f.setSize(4);
                readVerts(is, f.size(), nodeIdToPoints, verts);
                f.labelList::operator=(verts);
            }
        }
        else if (split[0] == "nsided")
        {
            label nShapes;
            is.read(nShapes);

            Pout<< indent << "nsided " << nShapes
                << " starting at line " << is.lineNumber()
                << " position " << is.stdStream().tellg() << endl;

            const label facei = faces.size();

            readIDs
            (
                is,
                read_elem_ids,
                nShapes,
                faceToElemIDs,
                elemIdToFaces
            );

            faces.setSize(facei+nShapes);

            for (label shapei = 0; shapei < nShapes; shapei++)
            {
                auto& f = faces[facei+shapei];
                label nVerts;
                is.read(nVerts);
                f.setSize(nVerts);
            }

            for (label shapei = 0; shapei < nShapes; shapei++)
            {
                auto& f = faces[facei+shapei];
                readVerts(is, f.size(), nodeIdToPoints, verts);
                f.labelList::operator=(verts);
            }
        }
        else
        {
            WarningInFunction << "Unhandled key " << string(split[0])
                << " from line " << line
                << " starting at line " << is.lineNumber()
                << " position " << is.stdStream().tellg() << endl;
        }
    }

    // EOF
    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool Foam::fileFormats::ensightMeshReader::readGeometry
(
    const scalar scaleFactor
)
{
    ensightReadFile is(geometryFile_);

    // Skip 'binary' tag
    is.readBinaryHeader();

    string header;
    is.read(header);
    Info<< "Ensight : " << header << endl;
    is.read(header);
    Info<< "Ensight : " << header << endl;


    bool read_node_ids = false;
    bool read_elem_ids = false;

    // Storage for all the parts
    // ~~~~~~~~~~~~~~~~~~~~~~~~~

    List<string> partNames;

    DynamicList<label> partIDs;     // per part the original Ensight part no
    PtrList<pointField> partPoints;
    PtrList<labelList> partNodeIDs;
    PtrList<Map<label>> partPointIDs;

    // Cells (cell-to-faces)
    // Note: only one internal mesh supported.
    PtrList<faceListList> partCells;
    // Element ID for cells
    PtrList<labelList> partCellIDs;
    PtrList<Map<label>> partCellElemIDs;
    // Boundary faces
    PtrList<faceList> partFaces;
    // Element IDs for faces
    PtrList<labelList> partFaceIDs;
    PtrList<Map<label>> partFaceElemIDs;


    // Parse all
    string line;
    while (is.good())
    {
        do
        {
            is.readKeyword(line);
        }
        while (line.empty() && is.good());
        const auto split = stringOps::splitSpace(line);

        if (split[0] == "extents")
        {
            point min;
            point max;
            is.read(min.x());
            is.read(max.x());
            is.read(min.y());
            is.read(max.y());
            is.read(min.z());
            is.read(max.z());
            Pout<< indent
                << "Read extents " << boundBox(min, max)
                << endl;
        }
        else if (split[0] == "node")
        {
            word id(split[1]);
            word op(split[2]);
            if (op == "given" || op == "ignore")
            {
                Pout<< indent << "Reading node ids" << endl;
                read_node_ids = true;
            }
        }
        else if (split[0] == "element")
        {
            word id(split[1]);
            word op(split[2]);
            if (op == "given" || op == "ignore")
            {
                Pout<< indent << "Reading element ids" << endl;
                read_elem_ids = true;
            }
        }
        else if (split[0] == "part")
        {
            bool finished = false;
            do
            {
                label partIndex;
                is.read(partIndex);
                partIDs.append(partIndex);

                // Make space
                partNames.setSize(partIDs.size());
                partPoints.append(new pointField(0));
                partNodeIDs.append(new labelList(0));
                partPointIDs.append(new Map<label>());

                partCells.append(new faceListList(0));
                partCellIDs.append(new labelList(0));
                partCellElemIDs.append(new Map<label>());

                partFaces.append(new faceList(0));
                partFaceIDs.append(new labelList(0));
                partFaceElemIDs.append(new Map<label>());

                is.read(partNames.last());

                Pout<< indent
                    << "Reading part " << partIndex
                    << " name " << partNames.last()
                    << " starting at line " << is.lineNumber()
                    << " position " << is.stdStream().tellg() << endl;

                Pout<< incrIndent;
                finished = readGoldPart
                (
                    is,
                    read_node_ids,
                    read_elem_ids,

                    partPoints.last(),
                    partNodeIDs.last(),
                    partPointIDs.last(),

                    // Cells (cell-to-faces)
                    partCells.last(),
                    partCellIDs.last(),
                    partCellElemIDs.last(),

                    // Faces
                    partFaces.last(),
                    partFaceIDs.last(),
                    partFaceElemIDs.last()
                );

                partPoints.last() *= scaleFactor;

                Pout<< indent
                    << "For part " << partIndex
                    << " read cells " << partCells.last().size()
                    << " faces " << partFaces.last().size()
                    << endl;

                Pout<< decrIndent;
            }
            while (!finished);

            break;
        }
    }


    // Merge all coordinates
    labelListList pointToMerged(partPoints.size());

    //- Use point merging - also merges points inside a part which might upset
    //- e.g. ami
    //{
    //    label nPoints = 0;
    //    forAll(partPoints, parti)
    //    {
    //        nPoints += partPoints[parti].size();
    //    }
    //
    //    points_.setSize(nPoints);
    //    nodeIds_.setSize(nPoints, -1);
    //    nPoints = 0;
    //    forAll(partPoints, parti)
    //    {
    //        const auto& pts = partPoints[parti];
    //
    //        SubList<point>(points_, pts.size(), nPoints) = pts;
    //        SubList<label>(nodeIds_, pts.size(), nPoints) =
    //            partNodeIDs[parti];
    //
    //        auto& map = pointToMerged[parti];
    //        map = nPoints + identity(pts.size());
    //
    //        nPoints += pts.size();
    //    }
    //
    //    if (mergeTol_ > 0)
    //    {
    //        const scalar mergeDist = mergeTol_*boundBox(points_, true).mag();
    //
    //        labelList sharedToMerged;
    //        const label nMerged = inplaceMergePoints
    //        (
    //            points_,
    //            mergeDist,
    //            false,
    //            sharedToMerged
    //        );
    //        Pout<< "Merged " << nMerged << " points out of " << nPoints
    //            << " using merge tolerance " << mergeTol_
    //            << ", distance " << mergeDist
    //            << endl;
    //
    //        forAll(partNodeIDs, parti)
    //        {
    //            inplaceRenumber(sharedToMerged, pointToMerged[parti]);
    //
    //            // Now pointToMerged contains the numbering from original,
    //            // partwise to global points
    //            UIndirectList<label>(nodeIds_, pointToMerged[parti]) =
    //                partNodeIDs[parti];
    //        }
    //    }
    //}

    // Merge all coordinates
    {
        boundBox extents;
        label nPoints = 0;
        forAll(partPoints, parti)
        {
            extents.add(partPoints[parti]);
            nPoints += partPoints[parti].size();
        }
        const scalar mergeDist = mergeTol_*extents.mag();

        Pout<< "Merging points closer than " << mergeDist
            << " calculated from bounding box " << extents
            << " and mergeTol " << mergeTol_
            << endl;

        forAll(partPoints, parti)
        {
            const auto& pPoints = partPoints[parti];
            auto& pPointMap = pointToMerged[parti];
            pPointMap.setSize(pPoints.size(), -1);

            Pout<< "Matching part " << parti
                << " name " << partNames[parti]
                << " points " << pPoints.size()
                << " to current merged points " << points_.size()
                << endl;

            if (points_.empty())
            {
                points_ = pPoints;
                pPointMap = identity(pPoints.size());
            }
            else
            {
                // Match to existing
                labelList from0To1;
                matchPoints
                (
                    pPoints,
                    points_,
                    scalarField(pPoints.size(), mergeDist),
                    false,
                    from0To1
                );
                label nAdded = 0;

                forAll(from0To1, partPointi)
                {
                    const label pointi = from0To1[partPointi];
                    if (pointi != -1)
                    {
                        pPointMap[partPointi] = pointi;
                    }
                    else
                    {
                        nAdded++;
                    }
                }

                const label nOldPoints = points_.size();
                points_.setSize(nOldPoints+nAdded);
                nAdded = 0;
                forAll(from0To1, partPointi)
                {
                    if (from0To1[partPointi] == -1)
                    {
                        const label newPointi = nOldPoints+nAdded++;
                        points_[newPointi] = pPoints[partPointi];
                        pPointMap[partPointi] = newPointi;
                    }
                }
            }
        }

        // Now we have complete points. Take over element IDs.
        nodeIds_.setSize(points_.size());
        forAll(partNodeIDs, parti)
        {
            const auto& pPointMap = pointToMerged[parti];
            UIndirectList<label>(nodeIds_, pPointMap) = partNodeIDs[parti];
        }

        Pout<< "Merged from " << nPoints << " down to " << points_.size()
            << " points" << endl;
    }


    // Merge all cells
    labelList cellOffsets(partCells.size()+1);
    cellOffsets[0] = 0;
    {
        label nCells = 0;
        forAll(partCells, parti)
        {
            nCells += partCells[parti].size();
            cellOffsets[parti+1] = nCells;
        }

        cellFaces_.setSize(nCells);
        cellTableId_.setSize(nCells);
        elementIds_.setSize(nCells, -1);

        forAll(partCells, parti)
        {
            // Copy cells into position
            const auto& cells = partCells[parti];

            SubList<faceList> cellSlice
            (
                cellFaces_,
                cells.size(),
                cellOffsets[parti]
            );
            cellSlice = cells;

            SubList<label>
            (
                cellTableId_,
                cells.size(),
                cellOffsets[parti]
            ) = parti;

            SubList<label>
            (
                elementIds_,
                cells.size(),
                cellOffsets[parti]
            ) = partCellIDs[parti];

            // Renumber faces
            const auto& pointMap = pointToMerged[parti];

            for (auto& cellFaces : cellSlice)
            {
                for (auto& f : cellFaces)
                {
                    inplaceRenumber(pointMap, f);
                }
            }
        }
    }


    // Add cells to separate zones
    forAll(partCells, parti)
    {
        cellTable_.setMaterial(parti, "fluid");
        cellTable_.setName(parti, "part" + Foam::name(partIDs[parti]));
    }


    // Build map from face to cell and index. Note: use exact match
    // - no circular equivalence
    // - but instead pass in ordered faces (lowest numbered vertex first)
    HashTable<cellFaceIdentifier, face, face::symmHasher> vertsToCell
    (
        2*cellOffsets.last()
    );

    // Insert cell's faces into hashtable
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    face rotatedFace;
    forAll(cellFaces_, celli)
    {
        const auto& cFaces = cellFaces_[celli];
        forAll(cFaces, cFacei)
        {
            const face& f = rotateFace(cFaces[cFacei], rotatedFace);

            const auto fFnd = vertsToCell.find(f);
            if (fFnd)
            {
                // Already inserted. Internal face.
                vertsToCell.erase(fFnd);
            }
            else
            {
                vertsToCell.insert(f, cellFaceIdentifier(celli, cFacei));
            }
        }
    }


    labelList patchToPart(partNames.size());
    label nPatches = 0;
    forAll(partFaces, parti)
    {
        if (partFaces[parti].size())
        {
            Pout<< "Using part " << parti
                << " name " << partNames[parti]
                << " as patch " << nPatches
                << endl;

            patchToPart[nPatches++] = parti;
        }
    }
    patchToPart.setSize(nPatches);
    boundaryIds_.setSize(nPatches);
    patchTypes_.setSize(nPatches, "wall");
    patchNames_.setSize(nPatches);
    forAll(patchNames_, patchi)
    {
        const label parti = patchToPart[patchi];

        Pout<< "Matching part " << parti
            << " name " << partNames[parti]
            << " faces " << partFaces[parti].size()
            //<< " points " << partPoints[parti].size()
            << " to merged faces " << vertsToCell.size()
            << ", merged points " << points_.size()
            << endl;

        patchNames_[patchi] = word::validate(partNames[parti]);

        const auto& pointMap = pointToMerged[parti];
        const auto& faces = partFaces[parti];

        auto& partCellAndFace = boundaryIds_[patchi];
        partCellAndFace.setSize(faces.size());

        label patchFacei = 0;
        forAll(faces, facei)
        {
            if (faces[facei].empty())
            {
                Pout<< "Ignoring empty face:" << facei << endl;
                continue;
            }

            // Rewrite into mesh-point ordering
            const face newF(pointMap, faces[facei]);
            // Lookup cell and face on cell using the vertices
            const auto cAndF = vertsToCell.find
            (
                rotateFace
                (
                    newF,
                    rotatedFace
                )
            );

            if (!cAndF)
            {
                //WarningInFunction
                //    << "Did not find face " << facei
                //    << " verts:" << faces[facei]
                //    << " renumbered:" << newF
                //    << " rotated:" << rotatedFace
                //    << " in part " << parti
                //    << endl;
            }
            else
            {
                partCellAndFace[patchFacei++] = cAndF();
                vertsToCell.erase(cAndF);
            }
        }
        partCellAndFace.setSize(patchFacei);
    }
    patchPhysicalTypes_.setSize(nPatches, "wall");
    patchStarts_.setSize(nPatches, 0);
    patchSizes_.setSize(nPatches, 0);

    if (vertsToCell.size())
    {
        // Unused internal or boundary faces
        boundaryIds_.append(List<cellFaceIdentifier>(0));
        {
            auto& cellAndFaces = boundaryIds_.last();
            cellAndFaces.setSize(vertsToCell.size());
            label i = 0;
            for (const auto& e : vertsToCell)
            {
                cellAndFaces[i++] = e;
            }
        }
        patchTypes_.append("empty");
        patchNames_.append("defaultFaces");
        patchPhysicalTypes_.append("empty");
        patchStarts_.append(0);
        patchSizes_.append(0);

        Pout<< "Introducing default patch " << patchNames_.size()-1
            << " name " << patchNames_.last() << endl;
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileFormats::ensightMeshReader::ensightMeshReader
(
    const fileName& geomFile,
    const objectRegistry& registry,
    const scalar mergeTol,
    const scalar scaleFactor,
    const bool setHandedness
)
:
    meshReader(geomFile, scaleFactor),
    mergeTol_(mergeTol),
    setHandedness_(setHandedness)
{}


// ************************************************************************* //
