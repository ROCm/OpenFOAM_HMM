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

Class
    Foam::foamVtkCells

Description
    The deep-copy description of an OpenFOAM volume mesh in data structures
    corresponding to an VTK UnstructuredGrid, including the possiblity of
    decomposing polyhedral cells into primitive cell types.

    Knowledge of the vtkUnstructuredGrid and the corresponding \c .vtu
    xml file-format aids in understanding this class.
    For flexibilty, support for the legacy vtk file-format is also provided.

    Primitive cell types are straighforward, polyhedral cells are represented
    by a face stream:
    \verbatim
        [nFaces, nFace0Pts, id1, id2, ..., nFace1Pts, id1, id2, ...]
    \endverbatim

    For the legacy format, the face stream is simply passed as vertex labels
    (connectivity).

    For the xml format, the face stream is saved separately:
    \verbatim
        "connectivity"
        == the unique vertex labels used by the cell (optionally sorted).

        "offsets":
        == offset + sizeof(connectivity)

        "faces":
        [nFaces, nFace0Pts, id1, id2, ..., nFace1Pts, id1, id2, ...]

        "faceoffsets":
        == faceoffsets + sizeof(faces)
    \endverbatim

    The storage of "connectivity" and "offsets" strongly resembles a
    CompactListList, but the "offsets" point to the end of the respective
    sub-lists.

SourceFiles
    foamVtkCells.C

\*---------------------------------------------------------------------------*/

#ifndef foamVtkCells_H
#define foamVtkCells_H

#include "foamVtkCore.H"
#include "DynamicList.H"
#include "SubList.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;

/*---------------------------------------------------------------------------*\
                      Class foamVtkCells Declaration
\*---------------------------------------------------------------------------*/

class foamVtkCells
:
    public fileFormats::foamVtkCore
{
public:

    //- Bookkeeping for polyhedral cell decomposition
    class decomp
    {
    private:
        friend foamVtkCells;

        // Private data

            //- Cell-centre labels for additional points of decomposed cells
            DynamicList<label> addPointCellLabels_;

            //- Label of original cell for decomposed cells
            DynamicList<label> superCells_;

            //- Point labels for subsetted meshes
            DynamicList<label> pointMap_;

            //- Track if decomposition was requested
            const bool requested_;


        // Private Member Functions

            //- Disallow default bitwise copy construct
            decomp(const decomp&) = delete;

            //- Disallow default bitwise assignment
            void operator=(const decomp&) = delete;

    public:

        // Constructors

            //- Construct null
            decomp(const bool decomposePoly = false);


        //- Destructor
        ~decomp();


        // Member Functions

        // Access

            //- Polyhedral decomposition requested
            inline bool requested() const;

            //- Polyhedral decomposition used
            inline bool used() const;

            //- Label of original cell for decomposed cells
            inline const labelList& superCells() const;

            //- Cell-centre labels for additional points of decomposed cells
            inline const labelList& addPointCellLabels() const;

            //- Point labels for subsetted meshes
            inline const labelList& pointMap() const;


        // Edit

            //- Clear
            void clear();
    };


private:

    // Private data

        //- Reference to underlying mesh or mesh sub-set
        const polyMesh& mesh_;

        //- Cell types (including added cells) in vtk numbering
        //  Range is 1-255
        List<uint8_t> cellTypes_;

        //- Vertices per cell (including added cells) in vtk ordering
        DynamicList<label> vertLabels_;

        //- Vertices per cell (including added cells) in vtk ordering
        DynamicList<label> vertOffset_;

        //- Face lists per polyhedral cell
        DynamicList<label> faceLabels_;

        //- Face label offsets
        DynamicList<label> faceOffset_;

        //- Bookkeeping for polyhedral cell decomposition
        decomp decompose_;

        //- Needs update
        bool needsUpdate_;



    // Private Member Functions

        //- Create the geometry
        void correct();

        //- Disallow default bitwise copy construct
        foamVtkCells(const foamVtkCells&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const foamVtkCells&) = delete;


public:

    // Constructors

        //- Construct from components.
        //  Optionally with polyhedral decomposition and/or lazy evaluation.
        //  A 'lazy' evaluation avoids fully creation within the constructor.
        foamVtkCells
        (
            const polyMesh&,
            const bool decomposePoly = false,
            const bool lazy = false
        );


    //- Destructor
    ~foamVtkCells();


    // Member Functions

        // Access

            //- Query the poly decompose flag.
            inline bool decomposeRequested() const;

            //- Values for "connectivity" (XML) or basis for "CELLS" (legacy)
            //  In the legacy format, the size (offset) must be prefixed.
            inline const labelList& vertLabels() const;

            //- Values for "offsets" (XML)
            //  or sizes to prefix for for "CELLS" (legacy)
            inline const labelList& vertOffsets() const;

            //- Values for "types" (XML) and "CELL_TYPES" (legacy)
            inline const List<uint8_t>& cellTypes() const;

            //- Values for "faces" (XML)
            inline const labelList& faceLabels() const;

            //- Values for "faceoffsets" (XML)
            inline const labelList& faceOffsets() const;

            //- Additional point addressing (from added point to original cell)
            inline const labelList& addPointCellLabels() const;

            //- Additional cells mapping (from added cell to original cell)
            inline const labelList& superCells() const;


            //- Number of field cells
            inline label nFieldCells() const;

            //- Number of field points
            label nFieldPoints() const;

            //- The field size for legacy "CELLS".
            //  In the legacy format, the size (offset) must be prefixed.
            label legacyCellPayLoad() const;


        //- Does the mapping need an update?
        bool needsUpdate() const;

        //- Mark as needing an update.
        //  May also free up unneeded data.
        //  Return false if it was already marked as expired.
        bool expire();

        //- Update the description (and decomposition) as required.
        //  Do nothing (and return false) if no update was required
        bool update();


        //- The const_iterator for foamVtkCells
        class const_iterator
        {
            friend class foamVtkCells;

        protected:

            // Protected Data

                //- Reference to parent list
                const foamVtkCells& parent_;

                //- Element index
                label index_;

                //- Begin of connectivity sub-list
                mutable label begVert_;

                //- Begin of faces sub-list
                mutable label begFace_;

                //- On-demand legacy pointer
                mutable autoPtr<SubList<label>> legacy_;


            // Constructors

                //- Construct begin/end iterator
                inline const_iterator
                (
                    const foamVtkCells&,
                    bool isEnd = false
                );

        public:

            // Member operators

                //- On-demand legacy cell labels (primitive or faces)
                inline const labelUList& legacyCell() const;

                //- Compare position
                inline bool operator!=(const const_iterator&) const;

                //- Pre-increment iterator
                inline const_iterator& operator++();
        };


        //- const_iterator set to the beginning
        inline const_iterator begin() const;

        //- const_iterator set to beyond the end
        inline const_iterator end() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "foamVtkCellsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
