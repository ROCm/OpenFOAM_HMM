/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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
#include "searchableSurfaces.H"
#include "polyMesh.H"
#include "IOdictionary.H"
#include "blockVertexList.H"
#include "blockEdgeList.H"
#include "blockFaceList.H"

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
    // Private Data

        //- Reference to mesh dictionary
        const IOdictionary& meshDict_;

        //- Switch for verbose output
        bool verboseOutput;

        //- Switch checking face consistency (defaults to true)
        bool checkFaceCorrespondence_;

        //- Optional searchable geometry to project face-points to
        searchableSurfaces geometry_;

        //- The scaling factor to convert to metres
        scalar scaleFactor_;

        //- The list of block vertices
        blockVertexList blockVertices_;

        //- The list of block vertex positions
        pointField vertices_;

        //- The list of curved edges
        blockEdgeList edges_;

        //- The list of curved faces
        blockFaceList faces_;

        //- The blocks themselves (the topology) as a polyMesh
        autoPtr<polyMesh> topologyPtr_;

        //- The sum of all cells in each block
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

        template<class Source>
        void checkPatchLabels
        (
            const Source& source,
            const word& patchName,
            const pointField& points,
            faceList& patchShapes
        ) const;

        void readPatches
        (
            const dictionary& meshDescription,
            faceListList& tmpBlocksPatches,
            wordList& patchNames,
            wordList& patchTypes,
            wordList& nbrPatchNames
        );

        void readBoundary
        (
            const dictionary& meshDescription,
            wordList& patchNames,
            faceListList& tmpBlocksPatches,
            PtrList<dictionary>& patchDicts
        );

        void createCellShapes(cellShapeList& tmpBlockCells);

        autoPtr<polyMesh> createTopology
        (
            const IOdictionary& dict,
            const word& regionName
        );

        void check(const polyMesh& bm, const dictionary& dict) const;

        //- Determine merge info and final number of cells/points
        //- based on point distances
        void calcGeometricalMerge();

        //- Determine merge info and final number of cells/points
        //- based on block topology
        void calcTopologicalMerge();

        faceList createPatchFaces(const polyPatch& patchTopologyFaces) const;

        void createPoints() const;
        void createCells() const;
        void createPatches() const;


        //- No copy construct
        blockMesh(const blockMesh&) = delete;

        //- No copy assignment
        void operator=(const blockMesh&) = delete;


public:

    // Data Types

        //- The block merging strategy
        enum mergeStrategy
        {
            TOPOLOGICAL,    //!< Merge using block topology
            GEOMETRIC       //!< Merge using point geometry
        };


    // Static Data Members

        ClassName("blockMesh");


    // Constructors

        //- Construct from IOdictionary for given region
        //  Default is topological merging.
        blockMesh
        (
            const IOdictionary& dict,
            const word& regionName,
            const mergeStrategy strategy = TOPOLOGICAL
        );


    //- Destructor
    ~blockMesh() = default;


    // Member Functions

        // Access

            //- Access to input dictionary
            const dictionary& meshDict() const
            {
                return meshDict_;
            }

            //- Optional searchable geometry to project face-points to
            const searchableSurfaces& geometry() const
            {
                return geometry_;
            }

            //- True if the blockMesh topology exists
            bool valid() const;

            //- Reference to point field defining the blockMesh
            //  these points have not been scaled by scaleFactor
            const pointField& vertices() const;

            //- Return the blockMesh topology as a polyMesh
            const polyMesh& topology() const;

            //- Return the curved edges
            const blockEdgeList& edges() const
            {
                return edges_;
            }

            //- Return the curved faces
            const blockFaceList& faces() const
            {
                return faces_;
            }

            //- The scaling factor used to convert to metres
            scalar scaleFactor() const;

            //- The points for the entire mesh
            //  these points have been scaled by scaleFactor
            const pointField& points() const;

            //- Return cell shapes list
            const cellShapeList& cells() const;

            //- Return the patch face lists
            const faceListList& patches() const;

            //- Get patch information from the topology mesh
            PtrList<dictionary> patchDicts() const;

            //- Return patch names
            wordList patchNames() const;

            //- Number of blocks with specified zones
            label numZonedBlocks() const;


        // Edit

            //- Enable/disable verbose information about the progress
            void verbose(const bool on=true);


        // Write

            //- Writes edges of blockMesh in OBJ format.
            void writeTopology(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
