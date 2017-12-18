/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenCFD Ltd.
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
    Foam::vtk::vtuCells

Description
    A deep-copy description of an OpenFOAM volume mesh in data structures
    suitable for VTK UnstructuredGrid, including the possibility of
    decomposing polyhedral cells into primitive cell types.

    Knowledge of the vtkUnstructuredGrid and the corresponding \c .vtu
    xml file-format aids in understanding this class.
    The class can be used for the VTK xml format, legacy format, as well as a
    VTK internal representation. The internal representation is somewhat
    related to the xml format, but not entirely.

SeeAlso
    Foam::vtk::vtuSizing

SourceFiles
    foamVtuCells.C
    foamVtuCellsI.H

\*---------------------------------------------------------------------------*/

#ifndef foamVtuCells_H
#define foamVtuCells_H

#include "foamVtkCore.H"
#include "foamVtkMeshMaps.H"
#include "foamVtuSizing.H"
#include "DynamicList.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;

namespace vtk
{
// Forward declaration of classes
    class outputOptions;

/*---------------------------------------------------------------------------*\
                     Class Foam::vtk::vtuCells Declaration
\*---------------------------------------------------------------------------*/

class vtuCells
:
    public vtuSizing
{
    // Private Member Data

      // Requested output types

        //- Output content type
        contentType output_;

        //- Bookkeeping for polyhedral cell decomposition
        bool decomposeRequest_;

      // Storage of output

        //- Cell types (including added cells) in vtk numbering
        //  Range is 1-255
        List<uint8_t> cellTypes_;

        //- Vertices per cell (including added cells) in vtk ordering
        List<label> vertLabels_;

        //- Connectivity (vertices) offset for the end of each cell
        List<label> vertOffset_;

        //- Face lists per polyhedral cell
        List<label> faceLabels_;

        //- Face label offsets
        List<label> faceOffset_;

        //- Bookkeeping for polyhedral cell decomposition
        foamVtkMeshMaps maps_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        vtuCells(const vtuCells&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const vtuCells&) = delete;


public:

    // Constructors

        //- Construct from components.
        //  Optionally with polyhedral decomposition.
        vtuCells
        (
            const enum contentType output = contentType::XML,
            const bool decompose = false
        );

        //- Construct from components and create the output information
        //  immediately
        vtuCells
        (
            const polyMesh& mesh,
            const enum contentType output = contentType::XML,
            const bool decompose = false
        );

        //- Construct from components.
        //  Optionally with polyhedral decomposition.
        vtuCells
        (
            const vtk::outputOptions outOpts,
            const bool decompose = false
        );

        //- Construct from components and create the output information
        //  immediately
        vtuCells
        (
            const polyMesh& mesh,
            const vtk::outputOptions outOpts,
            const bool decompose = false
        );


    //- Destructor
    ~vtuCells();


    // Member Functions

      // Access

        //- The output content type
        inline enum contentType content() const;

        //- Query the polyhedral decompose requested flag.
        inline bool decomposeRequested() const;

        //- True if no cellTypes are populated.
        inline bool empty() const;

        //- The size of populated cellTypes.
        inline label size() const;


      // Edit

        //- Reset all sizes to zero.
        void clear();

        //- Create the geometry using the previously requested output and
        //  decomposition types.
        void reset(const polyMesh& mesh);

        //- Respecify requested output and decomposition type prior to
        //  creating the geometry
        void reset
        (
            const polyMesh& mesh,
            const enum contentType output,
            const bool decompose
        );

        //- Renumber cell ids to account for subset meshes
        void renumberCells(const labelUList& mapping);

        //- Renumber point ids to account for subset meshes
        void renumberPoints(const labelUList& mapping);


      // Storage Access

        //- Values for "types" (XML) and "CELL_TYPES" (legacy)
        inline const List<uint8_t>& cellTypes() const;

        //- Values for "connectivity" (XML) or "CELLS" (legacy)
        inline const labelList& vertLabels() const;

        //- Values for "offsets" (XML only)
        inline const labelList& vertOffsets() const;

        //- Values for "faces" (XML only)
        inline const labelList& faceLabels() const;

        //- Values for "faceoffset" (XML only)
        inline const labelList& faceOffsets() const;

        //- Additional point addressing (from added point to original cell)
        inline const labelList& addPointCellLabels() const;

        //- Original cell ids for all cells (regular and decomposed).
        inline const labelList& cellMap() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace vtk
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "foamVtuCellsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
