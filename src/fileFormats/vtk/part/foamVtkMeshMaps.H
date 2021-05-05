/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2020 OpenCFD Ltd.
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
    Foam::foamVtkMeshMaps

Description
    Bookkeeping for mesh subsetting and/or polyhedral cell decomposition.
    Although the main use case is for handling vtk meshes, it is not specific
    to VTK alone.

    The cellMap is a local-to-global lookup for normal and decomposed cells.
    The pointMap is an optional local-to-global lookup for point ids.
    The additional ids is typically used to store the cell-centre labels
    for additional points of decomposed cells

SourceFiles
    foamVtkMeshMaps.C
    foamVtkMeshMapsI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_vtk_vtkMeshMaps_H
#define Foam_vtk_vtkMeshMaps_H

#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class foamVtkMeshMaps Declaration
\*---------------------------------------------------------------------------*/

class foamVtkMeshMaps
{
    // Private Data

        //- Original cell ids for all cells (regular and decomposed)
        DynamicList<label> cellMap_;

        //- Point labels for subsetted meshes
        DynamicList<label> pointMap_;

        //- Any additional (user) labels.
        //  Eg, cell-centre labels for additional points of decomposed cells
        DynamicList<label> additionalIds_;


public:

    // Constructors

        //- Default construct: zero-sized, no reserved size
        inline foamVtkMeshMaps();

        //- Construct with reserved size
        inline explicit foamVtkMeshMaps(const label size);


    // Member Functions

    // Access

        //- Original cell ids for all cells (regular and decomposed).
        //  A regular mesh comprising only primitive cell types, this will just
        //  be an identity list. However, for subsetted meshes and decomposed
        //  cells this becomes a useful means of mapping from the original mesh.
        inline const labelList& cellMap() const;

        //- Write access to original cell ids
        inline DynamicList<label>& cellMap();

        //- Point labels for subsetted meshes
        inline const labelList& pointMap() const;

        //- Write access to point labels for subsetted meshes
        inline DynamicList<label>& pointMap();

        //- Any additional (user) labels.
        //  Eg, cell-centre labels for additional points of decomposed cells
        inline const labelList& additionalIds() const;

        //- Write access to additional (user) labels.
        inline DynamicList<label>& additionalIds();


    // Edit

        //- Clear sizing
        inline void clear();

        //- Renumber cell ids (cellMap and additionalIds) to account for
        //  subset meshes
        void renumberCells(const labelUList& mapping);

        //- Renumber point ids (pointMap) to account for subset meshes
        void renumberPoints(const labelUList& mapping);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "foamVtkMeshMapsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
