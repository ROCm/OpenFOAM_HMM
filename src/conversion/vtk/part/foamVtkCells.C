/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenCFD Ltd.
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

#include "foamVtkCells.H"
#include "polyMesh.H"
#include "cellShape.H"
#include "cellModeller.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::foamVtkCells::correct()
{
    // Clear derived data
    // clearGeom();

    const cellModel& tet      = *(cellModeller::lookup("tet"));
    const cellModel& pyr      = *(cellModeller::lookup("pyr"));
    const cellModel& prism    = *(cellModeller::lookup("prism"));
    const cellModel& wedge    = *(cellModeller::lookup("wedge"));
    const cellModel& tetWedge = *(cellModeller::lookup("tetWedge"));
    const cellModel& hex      = *(cellModeller::lookup("hex"));

    const cellShapeList& cellShapes = mesh_.cellShapes();

    // face owner is needed to determine the face orientation
    const labelList& owner = mesh_.faceOwner();

    // Unique vertex labels per polyhedral
    HashSet<label> hashUniqId(2*256);

    // =======================
    // PASS 1: Determine sizes

    label nVertLabels = 0;
    label nFaceLabels = 0;
    label nAddPoints  = 0;
    label nAddCells   = 0;
    label nAddVerts   = 0;

    forAll(cellShapes, cellI)
    {
        const cellShape& shape = cellShapes[cellI];
        const cellModel& model = shape.model();

        if
        (
            model == tet
         || model == pyr
         || model == prism
         || model == hex
        )
        {
            // normal primitives
            nVertLabels += shape.size();
        }
        else if (model == tetWedge && decompose_.requested())
        {
            // Treat as squeezed prism (VTK_WEDGE)
            nVertLabels += 6;
        }
        else if (model == wedge && decompose_.requested())
        {
            // Treat as squeezed hex
            nVertLabels += 8;
        }
        else if (decompose_.requested())
        {
            // Polyhedral: Decompose into tets + pyramids.

            // Count vertices in first decomposed cell
            bool first = true;

            const cell& cFaces = mesh_.cells()[cellI];
            forAll(cFaces, cFaceI)
            {
                const face& f = mesh_.faces()[cFaces[cFaceI]];

                // Face decomposed into triangles and quads
                // Tri -> Tet, Quad -> Pyr
                label nTria = 0, nQuad = 0;
                f.nTrianglesQuads(mesh_.points(), nTria, nQuad);

                nAddCells  += nTria + nQuad;
                nAddVerts  += (nTria * 4) + (nQuad * 5);

                if (first)
                {
                    const label nvrt = (nQuad ? 5 : 4);
                    nAddCells--;
                    nAddVerts   -= nvrt;
                    nVertLabels += nvrt;

                    first = false;
                }
            }

            ++nAddPoints;
        }
        else
        {
            // Polyhedral: Not decomposed.

            const labelList& cFaces = mesh_.cells()[cellI];

            // establish unique node ids used (only needed for XML)
            hashUniqId.clear();

            // determing sizing for face stream
            // number of faces, size of each face, vertices per face
            // [nFaces, nFace0Pts, id1, id2, ..., nFace1Pts, id1, id2, ...]

            forAll(cFaces, cFaceI)
            {
                const face& f = mesh_.faces()[cFaces[cFaceI]];
                nFaceLabels += f.size();

                forAll(f, fp)
                {
                    hashUniqId.insert(f[fp]);
                }
            }

            nVertLabels += hashUniqId.size();
            nFaceLabels += 1 + cFaces.size();
        }
    }


    //
    // adjust/reserve sizes
    //

    // Cell types (including added cells) in vtk numbering
    cellTypes_.setSize(cellShapes.size() + nAddCells);

    // List of vertex labels in VTK ordering
    vertLabels_.setSize(nVertLabels + nAddVerts);

    vertOffset_.setSize(cellShapes.size() + nAddCells);

    faceLabels_.clear();
    faceOffset_.clear();
    if (nFaceLabels)
    {
        faceLabels_.setSize(nFaceLabels);

        // only need nCells (without nAddCells)
        // set to -1 (primitive)
        faceOffset_.setSize(cellShapes.size(), -1);
    }

    if (decompose_.requested())
    {
        decompose_.addPointCellLabels_.setSize(nAddPoints);
        decompose_.superCells_.setSize(nAddCells);
    }


    // ======================
    // PASS 2: Fill in arrays

    // Need this offset later, but only for decomposed polys
    const label offsetAddVerts = nVertLabels;

    // Reset counters
    nVertLabels = 0;
    nFaceLabels = 0;
    nAddPoints  = 0;
    nAddCells   = 0;
    nAddVerts   = 0;

    forAll(cellShapes, cellI)
    {
        const cellShape& shape = cellShapes[cellI];
        const cellModel& model = shape.model();

        if (model == tet)
        {
            cellTypes_[cellI] = foamVtkCore::VTK_TETRA;
            forAll(shape, i)
            {
                vertLabels_[nVertLabels++] = shape[i];
            }
            vertOffset_[cellI] = nVertLabels;
        }
        else if (model == pyr)
        {
            cellTypes_[cellI] = foamVtkCore::VTK_PYRAMID;
            forAll(shape, i)
            {
                vertLabels_[nVertLabels++] = shape[i];
            }
            vertOffset_[cellI] = nVertLabels;
        }
        else if (model == hex)
        {
            cellTypes_[cellI] = foamVtkCore::VTK_HEXAHEDRON;
            forAll(shape, i)
            {
                vertLabels_[nVertLabels++] = shape[i];
            }
            vertOffset_[cellI] = nVertLabels;
        }
        else if (model == prism)
        {
            cellTypes_[cellI] = foamVtkCore::VTK_WEDGE;

            // VTK_WEDGE triangles point outwards (swap 1<->2, 4<->5)
            vertLabels_[nVertLabels++] = shape[0];
            vertLabels_[nVertLabels++] = shape[2];
            vertLabels_[nVertLabels++] = shape[1];
            vertLabels_[nVertLabels++] = shape[3];
            vertLabels_[nVertLabels++] = shape[5];
            vertLabels_[nVertLabels++] = shape[4];

            vertOffset_[cellI] = nVertLabels;
        }
        else if (model == tetWedge && decompose_.requested())
        {
            // Treat as squeezed prism (VTK_WEDGE)
            cellTypes_[cellI] = foamVtkCore::VTK_WEDGE;

            vertLabels_[nVertLabels++] = shape[0];
            vertLabels_[nVertLabels++] = shape[2];
            vertLabels_[nVertLabels++] = shape[1];
            vertLabels_[nVertLabels++] = shape[3];
            vertLabels_[nVertLabels++] = shape[4];
            vertLabels_[nVertLabels++] = shape[3];

            vertOffset_[cellI] = nVertLabels;
        }
        else if (model == wedge && decompose_.requested())
        {
            // Treat as squeezed hex
            cellTypes_[cellI] = foamVtkCore::VTK_HEXAHEDRON;

            vertLabels_[nVertLabels++] = shape[0];
            vertLabels_[nVertLabels++] = shape[1];
            vertLabels_[nVertLabels++] = shape[2];
            vertLabels_[nVertLabels++] = shape[2];
            vertLabels_[nVertLabels++] = shape[3];
            vertLabels_[nVertLabels++] = shape[4];
            vertLabels_[nVertLabels++] = shape[5];
            vertLabels_[nVertLabels++] = shape[6];

            vertOffset_[cellI] = nVertLabels;
        }
        else if (decompose_.requested())
        {
            // Polyhedral cell - decompose into tet/pyr.

            // Ensure we have the correct orientation for the base of the
            // primitive cell shape.
            // If the cell is face owner, the orientation needs to be flipped
            // to avoid defining negative cells.
            // VTK doesn't seem to care, but we'll do it anyhow for safety.

            // The new vertex from the cell-centre
            const label newVertexLabel = mesh_.nPoints() + nAddPoints;

            // Mapping from additional point to cell
            decompose_.addPointCellLabels_[nAddPoints++] = cellI;

            // Whether to insert cell in place of original or not.
            bool first = true;

            const labelList& cFaces = mesh_.cells()[cellI];
            forAll(cFaces, cFaceI)
            {
                const face& f = mesh_.faces()[cFaces[cFaceI]];
                const bool isOwner = (owner[cFaces[cFaceI]] == cellI);

                // Count triangles/quads in decomposition
                label nTria = 0;
                label nQuad = 0;
                f.nTrianglesQuads(mesh_.points(), nTria, nQuad);

                // Do actual decomposition
                faceList faces3(nTria);
                faceList faces4(nQuad);
                nTria = 0, nQuad = 0;
                f.trianglesQuads(mesh_.points(), nTria, nQuad, faces3, faces4);

                forAll(faces4, fci)
                {
                    const face& quad = faces4[fci];

                    label celLoc;
                    label vrtLoc;

                    if (first)
                    {
                        celLoc = cellI;
                        vrtLoc = nVertLabels;
                        nVertLabels += 5;

                        vertOffset_[celLoc] = nVertLabels;
                        first = false;
                    }
                    else
                    {
                        celLoc = mesh_.nCells() + nAddCells;
                        vrtLoc = offsetAddVerts  + nAddVerts;
                        nAddVerts += 5;

                        vertOffset_[celLoc] = nAddVerts;
                        decompose_.superCells_[nAddCells++] = cellI;
                    }

                    cellTypes_[celLoc] = foamVtkCore::VTK_PYRAMID;

                    // See note above about the orientation.
                    if (isOwner)
                    {
                        vertLabels_[vrtLoc++] = quad[3];
                        vertLabels_[vrtLoc++] = quad[2];
                        vertLabels_[vrtLoc++] = quad[1];
                        vertLabels_[vrtLoc++] = quad[0];
                    }
                    else
                    {
                        vertLabels_[vrtLoc++] = quad[0];
                        vertLabels_[vrtLoc++] = quad[1];
                        vertLabels_[vrtLoc++] = quad[2];
                        vertLabels_[vrtLoc++] = quad[3];
                    }

                    vertLabels_[vrtLoc++] = newVertexLabel;
                }

                forAll(faces3, fci)
                {
                    const face& tria = faces3[fci];

                    label celLoc;
                    label vrtLoc;

                    if (first)
                    {
                        celLoc = cellI;
                        vrtLoc = nVertLabels;
                        nVertLabels += 4;

                        vertOffset_[celLoc] = nVertLabels;
                        first = false;
                    }
                    else
                    {
                        celLoc = mesh_.nCells() + nAddCells;
                        vrtLoc = offsetAddVerts + nAddVerts;
                        nAddVerts += 4;

                        vertOffset_[celLoc] = nAddVerts;
                        decompose_.superCells_[nAddCells++] = cellI;
                    }

                    cellTypes_[celLoc] = foamVtkCore::VTK_TETRA;

                    // See note above about the orientation.
                    if (isOwner)
                    {
                        vertLabels_[vrtLoc++] = tria[2];
                        vertLabels_[vrtLoc++] = tria[1];
                        vertLabels_[vrtLoc++] = tria[0];
                    }
                    else
                    {
                        vertLabels_[vrtLoc++] = tria[0];
                        vertLabels_[vrtLoc++] = tria[1];
                        vertLabels_[vrtLoc++] = tria[2];
                    }
                    vertLabels_[vrtLoc++] = newVertexLabel;
                }
            }
        }
        else
        {
            // Polyhedral cell - not decomposed

            hashUniqId.clear();  // unique node ids used (only needed for XML)

            // face-stream
            //   [nFaces, nFace0Pts, id1, id2, ..., nFace1Pts, id1, id2, ...]

            cellTypes_[cellI] = foamVtkCore::VTK_POLYHEDRON;
            const labelList& cFaces = mesh_.cells()[cellI];

            faceLabels_[nFaceLabels++] = cFaces.size();

            forAll(cFaces, cFaceI)
            {
                const face& f = mesh_.faces()[cFaces[cFaceI]];
                const bool isOwner = (owner[cFaces[cFaceI]] == cellI);

                forAll(f, fp)
                {
                    hashUniqId.insert(f[fp]);
                }

                // number of labels for this face
                faceLabels_[nFaceLabels++] = f.size();

                if (isOwner)
                {
                    forAll(f, fp)
                    {
                        faceLabels_[nFaceLabels++] = f[fp];
                    }
                }
                else
                {
                    // fairly immaterial if we reverse the list
                    // or use face::reverseFace()
                    forAllReverse(f, fp)
                    {
                        faceLabels_[nFaceLabels++] = f[fp];
                    }
                }
            }

            faceOffset_[cellI] = nFaceLabels;

            const labelList uniq = hashUniqId.sortedToc();
            forAll(uniq, i)
            {
                vertLabels_[nVertLabels++] = uniq[i];
            }

            vertOffset_[cellI] = nVertLabels;
        }
    }

    // ===========================================
    // PASS 3: Repair offsets for additional cells

//     Info<<"vertOffset: " << vertOffset_.size() << " VS. " << (mesh_.nCells()) << endl;
//     Info<<"nAddCells: "  << nAddCells << " VS. " << (mesh_.nCells()) << endl;

    if (nAddCells)
    {
        const label beg = mesh_.nCells();
        const label add = vertOffset_[beg-1];

        for (label i = beg; i < vertOffset_.size(); ++i)
        {
            vertOffset_[i] += add;
        }
    }

    // Some basic programming/sanity checks

    if ((nVertLabels + nAddVerts) != vertOffset_[mesh_.nCells()-1 + nAddCells])
    {
        WarningInFunction
            << "predicted offsets (" << nVertLabels << " + " << nAddVerts << ") != "
            << vertOffset_[mesh_.nCells()-1 + nAddCells]
            << endl;
    }

    if (offsetAddVerts != vertOffset_[mesh_.nCells()-1])
    {
        WarningInFunction
            << "predicted regular offset " << offsetAddVerts
            << " != " << vertOffset_[mesh_.nCells()]
            << endl;
    }

    // nFaceLabels = 0;
    // nAddPoints  = 0;
    // nAddCells   = 0;

    // Pout<<"vertLabels: " << vertLabels_.size() << " vs. " << (nVertLabels + nAddVerts) << endl;
    // Pout<<"faceLabels: " << faceLabels_.size() << " vs. " << nFaceLabels << endl;
#if 0
    if (decompose_.requested())
    {
        Pout<< "    Original cells:" << mesh_.nCells()
            << " points:" << mesh_.nPoints()
            << " Additional cells:" << decompose_.superCells_.size()
            << " additional points:" << decompose_.addPointCellLabels_.size()
            << nl << endl;
    }
#endif
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::foamVtkCells::decomp::decomp(const bool decomposePoly)
:
    addPointCellLabels_(),
    superCells_(),
    pointMap_(),
    requested_(decomposePoly)
{}


Foam::foamVtkCells::foamVtkCells
(
    const polyMesh& mesh,
    const bool decomposePoly,
    const bool lazy
)
:
    mesh_(mesh),
    cellTypes_(),
    vertLabels_(),
    vertOffset_(),
    faceLabels_(),
    faceOffset_(),
    decompose_(decomposePoly),
    needsUpdate_(true)
{
    if (!lazy)
    {
        correct();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::foamVtkCells::decomp::~decomp()
{}


Foam::foamVtkCells::~foamVtkCells()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::foamVtkCells::decomp::clear()
{
    superCells_.clear();
    addPointCellLabels_.clear();
    pointMap_.clear();
}


Foam::label Foam::foamVtkCells::nFieldPoints() const
{
    return mesh_.nPoints() + decompose_.addPointCellLabels_.size();
}


Foam::label Foam::foamVtkCells::legacyCellPayLoad() const
{
    label payLoad = cellTypes_.size();

    if (faceOffset_.size())
    {
        // also has polys with face streams

        label begVert = 0;
        label begFace = 0;

        forAll(faceOffset_, i)
        {
            label endFace = faceOffset_[i];
            label endVert = vertOffset_[i];

            if (endFace > 0)
            {
                // poly with face stream
                payLoad += endFace - begFace;

                begFace = endFace;
            }
            else
            {
                // primitive without face stream
                payLoad += endVert - begVert;
            }
            begVert = endVert;
        }
    }
    else if (vertOffset_.size())
    {
        // primitives only, trivial
        payLoad += vertOffset_[vertOffset_.size()-1];
    }

    return payLoad;
}


bool Foam::foamVtkCells::needsUpdate() const
{
    return needsUpdate_;
}


bool Foam::foamVtkCells::expire()
{
    // Clear any stored topologies

    // Clear derived data
    // clearGeom();

    // already marked as expired
    if (needsUpdate_)
    {
        return false;
    }

    needsUpdate_ = true;
    return true;
}


bool Foam::foamVtkCells::update()
{
    if (!needsUpdate_)
    {
        return false;
    }

    correct();

    needsUpdate_ = false;
    return true;
}


// ************************************************************************* //
