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

// Only used in this file
#include "foamVtuSizingTemplates.C"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::vtk::vtuSizing::presizeMaps(foamVtkMeshMaps& maps) const
{
    maps.cellMap().setSize(this->nFieldCells());
    maps.additionalIds().setSize(this->nAddPoints());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::vtk::vtuSizing::vtuSizing
(
    const polyMesh& mesh,
    const bool decompose
)
{
    clear();
    reset(mesh, decompose);
}


Foam::vtk::vtuSizing::vtuSizing()
{
    clear();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::vtk::vtuSizing::~vtuSizing()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::vtk::vtuSizing::reset
(
    const polyMesh& mesh,
    const bool decompose
)
{
    const cellModel& tet      = cellModel::ref(cellModel::TET);
    const cellModel& pyr      = cellModel::ref(cellModel::PYR);
    const cellModel& prism    = cellModel::ref(cellModel::PRISM);
    const cellModel& wedge    = cellModel::ref(cellModel::WEDGE);
    const cellModel& tetWedge = cellModel::ref(cellModel::TETWEDGE);
    const cellModel& hex      = cellModel::ref(cellModel::HEX);

    const cellShapeList& shapes = mesh.cellShapes();

    // Unique vertex labels per polyhedral
    HashSet<label> hashUniqId(2*256);

    decompose_ = decompose;
    nCells_    = mesh.nCells();
    nPoints_   = mesh.nPoints();
    nAddCells_ = 0;
    nAddVerts_ = 0;

    nCellsPoly_  = nCells_;
    nVertLabels_ = 0;
    nFaceLabels_ = 0;
    nVertPoly_   = 0;

    forAll(shapes, celli)
    {
        const cellShape& shape = shapes[celli];
        const cellModel& model = shape.model();

        if
        (
            model == tet
         || model == pyr
         || model == prism
         || model == hex
        )
        {
            // Normal primitive - not a poly
            --nCellsPoly_;
            nVertLabels_ += shape.size();
        }
        else if (model == tetWedge && decompose)
        {
            nVertLabels_ += 6;  // Treat as squeezed prism (VTK_WEDGE)
        }
        else if (model == wedge && decompose)
        {
            nVertLabels_ += 8;  // Treat as squeezed hex
        }
        else if (decompose)
        {
            // Polyhedral: Decompose into tets + pyramids.
            ++nAddPoints_;

            // Count vertices into first decomposed cell
            bool first = true;

            const cell& cFaces = mesh.cells()[celli];
            forAll(cFaces, cFaceI)
            {
                const face& f = mesh.faces()[cFaces[cFaceI]];

                // Face decomposed into triangles and quads
                // Tri -> Tet, Quad -> Pyr
                label nTria = 0, nQuad = 0;
                f.nTrianglesQuads(mesh.points(), nTria, nQuad);

                nAddCells_ += nTria + nQuad;
                nAddVerts_ += (nTria * 4) + (nQuad * 5);

                if (first)
                {
                    first = false;
                    --nAddCells_;

                    const label nvrt = (nQuad ? 5 : 4);
                    nAddVerts_   -= nvrt;
                    nVertLabels_ += nvrt;
                }
            }
        }
        else
        {
            // Polyhedral: Not decomposed

            const labelList& cFaces = mesh.cells()[celli];

            // Unique node ids used (XML/INTERNAL, not needed for LEGACY)
            hashUniqId.clear();

            // Face stream sizing:
            // number of faces, size of each face, vertices per face
            // [nFaces, nFace0Pts, id1, id2, ..., nFace1Pts, id1, id2, ...]

            forAll(cFaces, cFaceI)
            {
                const face& f = mesh.faces()[cFaces[cFaceI]];
                nFaceLabels_ += f.size();

                hashUniqId.insert(f);
            }

            // Legacy format only uses the face-stream.
            // - track what *NOT* to use for legacy
            nVertLabels_ += hashUniqId.size();
            nVertPoly_   += hashUniqId.size();

            nFaceLabels_ += 1 + cFaces.size();
        }
    }

    // decompose requested and needed
    decompose_ = (decompose && nCellsPoly_);
}


Foam::label Foam::vtk::vtuSizing::sizeOf
(
    const enum contentType output,
    const enum slotType slot
) const
{
    switch (output)
    {
        case contentType::LEGACY:
        {
            switch (slot)
            {
                case slotType::CELLS:
                    // legacy uses connectivity for primitives, but directly
                    // stores face streams into connectivity as well.
                    // size-prefix per cell
                    return
                    (
                        nVertLabels() + nAddVerts() - nVertPoly() // primitives
                      + nFaceLabels()     // face-stream (poly)
                      + nFieldCells()     // nFieldCells (size prefix)
                    );
                    break;

                default:
                    break;
            }
            break;
        }
        case contentType::XML:
        {
            switch (slot)
            {
                case slotType::CELLS:
                    return (nVertLabels() + nAddVerts());
                    break;

                case slotType::CELLS_OFFSETS:
                    return nFieldCells();
                    break;

                case slotType::FACES:
                    return nFaceLabels();
                    break;

                case slotType::FACES_OFFSETS:
                    return nFaceLabels() ? nFieldCells() : 0;
                    break;
            }
            break;
        }
        case contentType::INTERNAL:
        {
            switch (slot)
            {
                case slotType::CELLS:
                    // size-prefix per cell
                    return (nVertLabels() + nAddVerts() + nFieldCells());
                    break;

                case slotType::CELLS_OFFSETS:
                    return nFieldCells();
                    break;

                case slotType::FACES:
                    return nFaceLabels();
                    break;

                case slotType::FACES_OFFSETS:
                    return nFaceLabels() ? nFieldCells() : 0;
                    break;
            }
            break;
        }
    }

    return 0;
}


void Foam::vtk::vtuSizing::populateLegacy
(
    const polyMesh& mesh,
    UList<uint8_t>& cellTypes,
    labelUList& vertLabels,
    foamVtkMeshMaps& maps
) const
{
    // Leave as zero-sized so that populateArrays doesn't fill it.
    List<label> unused;

    presizeMaps(maps);

    populateArrays
    (
        mesh,
        *this,
        cellTypes,
        vertLabels,
        unused, // offsets
        unused, // faces
        unused, // facesOffsets
        contentType::LEGACY,
        maps.cellMap(),
        maps.additionalIds()
    );
}


void Foam::vtk::vtuSizing::populateXml
(
    const polyMesh& mesh,
    UList<uint8_t>& cellTypes,
    labelUList& connectivity,
    labelUList& offsets,
    labelUList& faces,
    labelUList& facesOffsets,
    foamVtkMeshMaps& maps
) const
{
    presizeMaps(maps);

    populateArrays
    (
        mesh,
        *this,
        cellTypes,
        connectivity,
        offsets,
        faces,
        facesOffsets,
        contentType::XML,
        maps.cellMap(),
        maps.additionalIds()
    );
}


void Foam::vtk::vtuSizing::populateInternal
(
    const polyMesh& mesh,
    UList<uint8_t>& cellTypes,
    UList<int>& connectivity,
    UList<int>& offsets,
    UList<int>& faces,
    UList<int>& facesOffsets,
    foamVtkMeshMaps& maps
) const
{
    presizeMaps(maps);

    populateArrays
    (
        mesh,
        *this,
        cellTypes,
        connectivity,
        offsets,
        faces,
        facesOffsets,
        contentType::INTERNAL,
        maps.cellMap(),
        maps.additionalIds()
    );
}


void Foam::vtk::vtuSizing::populateInternal
(
    const polyMesh& mesh,
    UList<uint8_t>& cellTypes,
    UList<long>& connectivity,
    UList<long>& offsets,
    UList<long>& faces,
    UList<long>& facesOffsets,
    foamVtkMeshMaps& maps
) const
{
    presizeMaps(maps);

    populateArrays
    (
        mesh,
        *this,
        cellTypes,
        connectivity,
        offsets,
        faces,
        facesOffsets,
        contentType::INTERNAL,
        maps.cellMap(),
        maps.additionalIds()
    );
}


void Foam::vtk::vtuSizing::populateInternal
(
    const polyMesh& mesh,
    UList<uint8_t>& cellTypes,
    UList<long long>& connectivity,
    UList<long long>& offsets,
    UList<long long>& faces,
    UList<long long>& facesOffsets,
    foamVtkMeshMaps& maps
) const
{
    presizeMaps(maps);

    populateArrays
    (
        mesh,
        *this,
        cellTypes,
        connectivity,
        offsets,
        faces,
        facesOffsets,
        contentType::INTERNAL,
        maps.cellMap(),
        maps.additionalIds()
    );
}


void Foam::vtk::vtuSizing::populateInternal
(
    const polyMesh& mesh,
    UList<uint8_t>& cellTypes,
    UList<int>& connectivity,
    UList<int>& offsets,
    UList<int>& faces,
    UList<int>& facesOffsets,
    labelUList& cellMap,
    labelUList& addPointsIds
) const
{
    populateArrays
    (
        mesh,
        *this,
        cellTypes,
        connectivity,
        offsets,
        faces,
        facesOffsets,
        contentType::INTERNAL,
        cellMap,
        addPointsIds
    );
}


void Foam::vtk::vtuSizing::populateInternal
(
    const polyMesh& mesh,
    UList<uint8_t>& cellTypes,
    UList<long>& connectivity,
    UList<long>& offsets,
    UList<long>& faces,
    UList<long>& facesOffsets,
    labelUList& cellMap,
    labelUList& addPointsIds
) const
{
    populateArrays
    (
        mesh,
        *this,
        cellTypes,
        connectivity,
        offsets,
        faces,
        facesOffsets,
        contentType::INTERNAL,
        cellMap,
        addPointsIds
    );
}


void Foam::vtk::vtuSizing::populateInternal
(
    const polyMesh& mesh,
    UList<uint8_t>& cellTypes,
    UList<long long>& connectivity,
    UList<long long>& offsets,
    UList<long long>& faces,
    UList<long long>& facesOffsets,
    labelUList& cellMap,
    labelUList& addPointsIds
) const
{
    populateArrays
    (
        mesh,
        *this,
        cellTypes,
        connectivity,
        offsets,
        faces,
        facesOffsets,
        contentType::INTERNAL,
        cellMap,
        addPointsIds
    );
}


void Foam::vtk::vtuSizing::clear()
{
    decompose_   = false;
    nCells_      = 0;
    nPoints_     = 0;
    nVertLabels_ = 0;

    nFaceLabels_ = 0;
    nCellsPoly_  = 0;
    nVertPoly_   = 0;

    nAddCells_   = 0;
    nAddPoints_  = 0;
    nAddVerts_   = 0;
}


void Foam::vtk::vtuSizing::info(Ostream& os) const
{
    os  << "nFieldCells:" << nFieldCells();
    if (nAddCells_)
    {
        os  << " (" << nCells_
            << "+" << nAddCells_ << ")";
    }
    else
    {
        os  << " (poly:" << nCellsPoly_ << ")";
    }

    os  << " nFieldPoints:" << nFieldPoints();
    if (nAddPoints_)
    {
        os  << " (" << nPoints_ << "+" << nAddPoints_ << ")";
    }

    os  << " nVertLabels:" << (nVertLabels_ + nAddVerts_);
    if (nAddVerts_)
    {
        os  << " (" << nVertLabels_ << "+" << nAddVerts_ << ")";
    }
    else if (nVertPoly_)
    {
        os  << " (poly:" << nVertPoly_ << ")";
    }

    os << " nFaceLabels:" << nFaceLabels_;
    os << " legacy-count:" << sizeLegacy();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool Foam::vtk::vtuSizing::operator==(const vtuSizing& rhs) const
{
    return
    (
        decompose()   == rhs.decompose()
     && nCells()      == rhs.nCells()
     && nPoints()     == rhs.nPoints()
     && nVertLabels() == rhs.nVertLabels()
     && nFaceLabels() == rhs.nFaceLabels()
     && nCellsPoly()  == rhs.nCellsPoly()
     && nVertPoly()   == rhs.nVertPoly()
     && nAddCells()   == rhs.nAddCells()
     && nAddPoints()  == rhs.nAddPoints()
     && nAddVerts()   == rhs.nAddVerts()
    );
}


bool Foam::vtk::vtuSizing::operator!=(const vtuSizing& rhs) const
{
    return !operator==(rhs);
}


// ************************************************************************* //
