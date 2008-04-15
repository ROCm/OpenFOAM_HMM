/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

Class
    Foam::mapDistributePolyMesh

Description
    Class containing mesh-to-mesh mapping information after a mesh distribution
    where we send parts of meshes (using subsetting) to other processors 
    and receive and reconstruct mesh.

    We store mapping from the bits-to-send to the complete starting mesh
    (subXXXMap) and from the received bits to their location in the new
    mesh (constructXXXMap).

SourceFiles
    mapDistributePolyMesh.C

\*---------------------------------------------------------------------------*/

#ifndef mapDistributePolyMesh_H
#define mapDistributePolyMesh_H

#include "labelList.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                           Class mapDistributePolyMesh Declaration
\*---------------------------------------------------------------------------*/

class mapDistributePolyMesh
{
    // Private data

        const polyMesh& mesh_;

        //- Number of old live points
        const label nOldPoints_;

        //- Number of old live faces
        const label nOldFaces_;

        //- Number of old live cells
        const label nOldCells_;

        //- List of the old patch sizes
        labelList oldPatchSizes_;

        //- List of the old patch start labels
        const labelList oldPatchStarts_;

        //- List of numbers of mesh points per old patch
        const labelList oldPatchNMeshPoints_;

        //- Maps from subsetted mesh back to original mesh
        const labelListList subPointMap_;
        const labelListList subFaceMap_;
        const labelListList subCellMap_;
        const labelListList subPatchMap_;

        //- Maps from subsetted mesh to new reconstructed mesh
        const labelListList constructPointMap_;
        const labelListList constructFaceMap_;
        const labelListList constructCellMap_;
        const labelListList constructPatchMap_;

    // Private Member Functions

        void calcPatchSizes();


public:

    // Constructors

        //- Construct from components
        mapDistributePolyMesh
        (
            const polyMesh& mesh,

            // mesh before changes
            const label nOldPoints,
            const label nOldFaces,
            const label nOldCells,
            const labelList& oldPatchStarts,
            const labelList& oldPatchNMeshPoints,

            // how to subset pieces of mesh to send across
            const labelListList& subPointMap,
            const labelListList& subFaceMap,
            const labelListList& subCellMap,
            const labelListList& subPatchMap,

            // how to reconstruct received mesh
            const labelListList& constructPointMap,
            const labelListList& constructFaceMap,
            const labelListList& constructCellMap,
            const labelListList& constructPatchMap
        );

        //- (optionally destructively) construct from components
        mapDistributePolyMesh
        (
            const polyMesh& mesh,
            const label nOldPoints,
            const label nOldFaces,
            const label nOldCells,
            labelList& oldPatchStarts,
            labelList& oldPatchNMeshPoints,

            labelListList& subPointMap,
            labelListList& subFaceMap,
            labelListList& subCellMap,
            labelListList& subPatchMap,
            labelListList& constructPointMap,
            labelListList& constructFaceMap,
            labelListList& constructCellMap,
            labelListList& constructPatchMap,
            const bool reUse                // clone or reuse
        );


    // Member Functions

        // Access

            const polyMesh& mesh() const
            {
                return mesh_;
            }

            //- Number of points in mesh before distribution
            label nOldPoints() const
            {
                return nOldPoints_;
            }

            //- Number of faces in mesh before distribution
            label nOldFaces() const
            {
                return nOldFaces_;
            }

            //- Number of cells in mesh before distribution
            label nOldCells() const
            {
                return nOldCells_;
            }

            //- List of the old patch sizes
            const labelList& oldPatchSizes() const
            {
                return oldPatchSizes_;
            }

            //- List of the old patch start labels
            const labelList& oldPatchStarts() const
            {
                return oldPatchStarts_;
            }

            //- List of numbers of mesh points per old patch
            const labelList& oldPatchNMeshPoints() const
            {
                return oldPatchNMeshPoints_;
            }

            //- From subsetted mesh back to original mesh
            //  (these are like the maps in fvMeshSubset)
            const labelListList& subPointMap() const
            {
                return subPointMap_;
            }
            //- From subsetted mesh back to original mesh
            const labelListList& subFaceMap() const
            {
                return subFaceMap_;
            }
            //- From subsetted mesh back to original mesh
            const labelListList& subCellMap() const
            {
                return subCellMap_;
            }
            //- From subsetted mesh back to original mesh
            const labelListList& subPatchMap() const
            {
                return subPatchMap_;
            }

            //- From subsetted mesh to new reconstructed mesh
            //  (like the pointProcAddressing in reconstructPar)
            const labelListList& constructPointMap() const
            {
                return constructPointMap_;
            }
            //- From subsetted mesh to new reconstructed mesh
            const labelListList& constructFaceMap() const
            {
                return constructFaceMap_;
            }
            //- From subsetted mesh to new reconstructed mesh
            const labelListList& constructCellMap() const
            {
                return constructCellMap_;
            }
            //- From subsetted mesh to new reconstructed mesh
            const labelListList& constructPatchMap() const
            {
                return constructPatchMap_;
            }

        // Edit

            //- Helper function: distribute list data
            template<class T>
            static void distribute
            (
                const label newSize,
                const labelListList& subMap,        // from subset to original
                const labelListList& constructMap,  // from subset to new
                List<T>&
            );

            //- distribute list of point data
            template<class T>
            void distributePointData(List<T>& lst) const
            {
                distribute
                (
                    mesh_.nPoints(),
                    subPointMap(),
                    constructPointMap(),
                    lst
                );
            }
            //- distribute list of face data
            template<class T>
            void distributeFaceData(List<T>& lst) const
            {
                distribute
                (
                    mesh_.nFaces(),
                    subFaceMap(),
                    constructFaceMap(),
                    lst
                );
            }
            //- distribute list of cell data
            template<class T>
            void distributeCellData(List<T>& lst) const
            {
                distribute
                (
                    mesh_.nCells(),
                    subCellMap(),
                    constructCellMap(),
                    lst
                );
            }
            //- distribute list of patch data
            template<class T>
            void distributePatchData(List<T>& lst) const
            {
                distribute
                (
                    mesh_.boundaryMesh().size(),
                    subPatchMap(),
                    constructPatchMap(),
                    lst
                );
            }


            //- distribute list of point/face/cell/patch indices.
            //  (Converts to boolList, distributes boolList and reconstructs)
            void distributePointIndices(labelList& pointIDs) const;
            void distributeFaceIndices(labelList& faceIDs) const;
            void distributeCellIndices(labelList& cellIDs) const;
            void distributePatchIndices(labelList& patchIDs) const;


            //- Correct for topo change.
            void updateMesh(const mapPolyMesh&)
            {
                notImplemented
                (
                    "mapDistributePolyMesh::updateMesh(const mapPolyMesh&)"
                );
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "mapDistributePolyMeshTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
