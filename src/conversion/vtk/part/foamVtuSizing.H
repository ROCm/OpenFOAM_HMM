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

Class
    Foam::vtk::vtuSizing

Description
    Sizing descriptions and routines for transcribing an OpenFOAM volume mesh
    into a VTK unstructured grid, with possible decomposition of polyhedral
    cells into primitive cell types.

    This class is intended to populate externally allocated arrays with content
    that is compatible with what VTK expects. This approach allows an improved
    separation of the OpenFOAM mesh description and the storage, and allows
    support of alternative storage containers (eg, std::vector, vtkDataArray).
    The ideal goal would be a zero-copy mechanism, but this does not work for
    several reasons:
    \par
    - OpenFOAM and VTK have different point ordering for prism
    - polyhedral decomposition
    - face-stream are required for VTK
    - VTK internal storage includes list size as part of the data
    - VTK includes storage may be a different base size (eg, long long)
      compared to the OpenFOAM label.

    \par Data Entries (slots)

    These are the storage entries normally associate with each output-type:
    \table
        legacy output
        \c types    | vtk cell type (1-255)
        \c cells    | nLabels and unique vertex labels used by the cell, or
                    | [nLabels nFaces, nFace0Pts, id1, id2, ..., nFace1Pts, id1, id2, ...]
    \endtable

    \table
        xml output
        \c types        | vtk cell type (1-255)
        \c connectivity | unique vertex labels used by the cell
        \c offsets      | end offset for each of \c connectivity
        \c faces        | face stream for polyhedral cells
                        | [nFaces, nFace0Pts, id1, id2, ..., nFace1Pts, id1, id2, ...]
        \c faceoffsets  | end offset for each of \c faces, with -1 for primitive cells
    \endtable

    \table
        internal storage
        \c types        | vtk cell type (1-255)
        \c connectivity | nLabels and unique vertex labels used by the cell
        \c location     | begin location for each of \c connectivity
        \c faces        | face stream for polyhedral cells
                        | [nFaces, nFace0Pts, id1, id2, ..., nFace1Pts, id1, id2, ...]
        \c facelocation | begin location for each of \c faces, with -1 for primitive cells
    \endtable

    The VTK storage concept for "connectivity" and "faces" somewhat resemble
    a CompactListList.

SourceFiles
    foamVtuSizing.C
    foamVtuSizingI.H

\*---------------------------------------------------------------------------*/

#ifndef foamVtuSizing_H
#define foamVtuSizing_H

#include "label.H"
#include "labelList.H"
#include "foamVtkMeshMaps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;

namespace vtk
{

/*---------------------------------------------------------------------------*\
                    Class Foam::vtk::vtuSizing Declaration
\*---------------------------------------------------------------------------*/

class vtuSizing
{
public:

    // Public data

    //- Types of content that the storage may represent
    enum contentType
    {
        LEGACY,     //!< Legacy VTK content
        XML,        //!< XML (VTU) content
        INTERNAL    //!< Internal vtkUnstructuredGrid content
    };

    //- The possible storage 'slots' that can be used
    enum slotType
    {
        CELLS,         //!< Cell connectivity (ALL)
        CELLS_OFFSETS, //!< End-offsets (XML) or locations (INTERNAL) for cells
        FACES,         //!< Face-stream (XML, INTERNAL)
        FACES_OFFSETS  //!< End-offsets (XML) or locations (INTERNAL) for faces
    };


private:

    // Private Member Data

        //- Polyhedral decomposition requested
        bool decompose_;

        //- Number of cells in the mesh
        label nCells_;

        //- Number of points in the mesh
        label nPoints_;

        //- Number of vertex labels to represent the mesh
        label nVertLabels_;

      // Polyhedrals

        //- Number of polyhedral face labels for the mesh
        label nFaceLabels_;

        //- Number of polyhedral cells (informational)
        label nCellsPoly_;

        //- Number of vertex labels used by polyhedrals
        label nVertPoly_;

      // Decomposed polyhedrals

        //- Number of additional (decomposed) cells for the mesh
        label nAddCells_;

        //- Number of additional (decomposed) points for the mesh
        label nAddPoints_;

        //- Number of additional (decomposed) vertices for the mesh
        label nAddVerts_;


    // Private Member Functions

        //- set-size for cellMap and additionalIds
        void presizeMaps(foamVtkMeshMaps& maps) const;

        //- Populate lists for internal VTK representation
        template<class LabelType, class LabelType2>
        static void populateArrays
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
        );

        // Allow default bitwise copy/assignment

public:

    // Constructors

        //- Construct null.
        vtuSizing();

        //- Construct sizes by analyzing the mesh,
        //  optionally with polyhedral decomposition.
        vtuSizing(const polyMesh& mesh, const bool decompose=false);


    //- Destructor
    ~vtuSizing();


    // Member Functions

      // Edit

        //- Construct sizes by analyzing the mesh,
        //  optionally with polyhedral decomposition.
        void reset(const polyMesh& mesh, const bool decompose=false);

        //- Reset all sizes to zero.
        void clear();


      // Access

        //- Query the decompose flag
        inline bool decompose() const;

        //- Number of cells for the mesh
        inline label nCells() const;

        //- Number of points for the mesh
        inline label nPoints() const;

        //- Number of vertex labels for the mesh
        inline label nVertLabels() const;

        //- Number of polyhedral face labels for the mesh
        inline label nFaceLabels() const;

        //- Number of polyhedral cells for the mesh
        inline label nCellsPoly() const;

        //- Number of vertex labels for polyhedral cells of the mesh
        inline label nVertPoly() const;

        //- Number of additional (decomposed) cells for the mesh
        inline label nAddCells() const;

        //- Number of additional (decomposed) points for the mesh
        inline label nAddPoints() const;

        //- Number of additional (decomposed) vertices for the mesh
        inline label nAddVerts() const;


        //- Number of field cells = nCells + nAddCells
        inline label nFieldCells() const;

        //- Number of field points = nPoints + nAddPoints
        inline label nFieldPoints() const;


      // Derived sizes

        //- Return the required size for the storage slot
        label sizeOf
        (
            const enum contentType output,
            const enum slotType slot
        ) const;


        //- The calculated size for legacy storage
        inline label sizeLegacy() const;

        //- The calculated size for legacy storage of the specified slot
        inline label sizeLegacy(const enum slotType slot) const;

        //- The calculated size for xml storage of the specified slot
        inline label sizeXml(const enum slotType slot) const;

        //- The calculated size for vtk-internal storage of the specified slot
        inline label sizeInternal(const enum slotType slot) const;


      // Utilty routines

        //- Populate lists for Legacy output
        void populateLegacy
        (
            const polyMesh& mesh,
            UList<uint8_t>& cellTypes,
            labelUList& connectivity,
            foamVtkMeshMaps& maps
        ) const;

        //- Populate lists for XML output
        void populateXml
        (
            const polyMesh& mesh,
            UList<uint8_t>& cellTypes,
            labelUList& connectivity,
            labelUList& offsets,
            labelUList& faces,
            labelUList& facesOffsets,
            foamVtkMeshMaps& maps
        ) const;

        //- Populate lists for Internal VTK format
        void populateInternal
        (
            const polyMesh& mesh,
            UList<uint8_t>& cellTypes,
            UList<int>& connectivity,
            UList<int>& offsets,
            UList<int>& faces,
            UList<int>& facesOffsets,
            foamVtkMeshMaps& maps
        ) const;

        //- Populate lists for Internal VTK format
        void populateInternal
        (
            const polyMesh& mesh,
            UList<uint8_t>& cellTypes,
            UList<long>& connectivity,
            UList<long>& offsets,
            UList<long>& faces,
            UList<long>& facesOffsets,
            foamVtkMeshMaps& maps
        ) const;

        //- Populate lists for Internal VTK format
        void populateInternal
        (
            const polyMesh& mesh,
            UList<uint8_t>& cellTypes,
            UList<long long>& connectivity,
            UList<long long>& offsets,
            UList<long long>& faces,
            UList<long long>& facesOffsets,
            foamVtkMeshMaps& maps
        ) const;

        //- Populate lists for Internal VTK format
        void populateInternal
        (
            const polyMesh& mesh,
            UList<uint8_t>& cellTypes,
            UList<int>& connectivity,
            UList<int>& offsets,
            UList<int>& faces,
            UList<int>& facesOffsets,
            labelUList& cellMap,
            labelUList& addPointsIds
        ) const;

        //- Populate lists for Internal VTK format
        void populateInternal
        (
            const polyMesh& mesh,
            UList<uint8_t>& cellTypes,
            UList<long>& connectivity,
            UList<long>& offsets,
            UList<long>& faces,
            UList<long>& facesOffsets,
            labelUList& cellMap,
            labelUList& addPointsIds
        ) const;

        //- Populate lists for Internal VTK format
        void populateInternal
        (
            const polyMesh& mesh,
            UList<uint8_t>& cellTypes,
            UList<long long>& connectivity,
            UList<long long>& offsets,
            UList<long long>& faces,
            UList<long long>& facesOffsets,
            labelUList& cellMap,
            labelUList& addPointsIds
        ) const;


      // Write

        //- Report some information
        void info(Ostream& os) const;


      // Member Operators

        //- Test equality
        bool operator==(const vtuSizing& rhs) const;

        //- Test inequality
        bool operator!=(const vtuSizing& rhs) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "foamVtuSizingI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
