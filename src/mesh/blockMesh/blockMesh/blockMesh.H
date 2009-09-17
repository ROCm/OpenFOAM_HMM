/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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
    Foam::blockMesh

Description
    A multi-block mesh generator

Note
    The vertices, cells and patches for filling the blocks are demand-driven.

SourceFiles
    blockMesh.C
    blockMeshCheck.C
    blockMeshCreate.C
    blockMeshMerge.C
    blockMeshTopology.C

\*---------------------------------------------------------------------------*/

#ifndef blockMesh_H
#define blockMesh_H

#include "blockList.H"
#include "polyMesh.H"
#include "IOdictionary.H"
#include "curvedEdgeList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class blockMesh Declaration
\*---------------------------------------------------------------------------*/

class blockMesh
:
    public blockList
{
    // Private data

        //- Point field defining the block mesh (corners)
        pointField blockPointField_;

        //- The list of curved edges
        curvedEdgeList edges_;

        //- The scaling factor to convert to meters
        scalar scaleFactor_;

        //- The blocks themselves (the topology) as a polyMesh
        polyMesh* topologyPtr_;

        label nPoints_;

        //- The sum of all cells in each block
        label nCells_;

        //- The point offset added to each block
        labelList blockOffsets_;

        //- The merge points information
        labelList mergeList_;

        mutable pointField points_;

        mutable cellShapeList cells_;

        mutable faceListList patches_;


    // Private Member Functions

        bool blockLabelsOK
        (
            const label blockLabel,
            const pointField& points,
            const cellShape& blockShape
        ) const;

        bool patchLabelsOK
        (
            const label patchLabel,
            const pointField& points,
            const faceList& patchShapes
        ) const;

        polyMesh* createTopology(IOdictionary&);
        void checkBlockMesh(const polyMesh&) const;

        //- Determine the merge info and the final number of cells/points
        void calcMergeInfo();

        faceList createPatchFaces(const polyPatch& patchTopologyFaces) const;

        void createPoints() const;
        void createCells() const;
        void createPatches() const;

        //- as copy (not implemented)
        blockMesh(const blockMesh&);


public:

    // Constructors

        //- Construct from IOdictionary
        blockMesh(IOdictionary&);


    // Destructor

        ~blockMesh();


    // Member Functions

        // Access

            const polyMesh& topology() const;

            const curvedEdgeList& edges() const
            {
                return edges_;
            }

            const pointField& points() const;

            const cellShapeList& cells() const;

            const faceListList& patches() const;

            wordList patchNames() const;

            wordList patchTypes() const;

            wordList patchPhysicalTypes() const;

            //- Number of blocks with specified zones
            label numZonedBlocks() const;


        // Write

            //- Writes edges of blockMesh in OBJ format.
            void writeTopology(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

