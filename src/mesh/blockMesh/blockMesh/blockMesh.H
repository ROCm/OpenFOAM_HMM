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

    Dictionary controls
    \table
        Property    | Description                           | Required | Default
        scale       | Point scaling                         | no  | 1.0
        vertices    |                                       | yes |
        blocks      |                                       | yes |
        edges       |                                       | no  |
        faces       |                                       | no  |
        boundary    | Boundary definition                   | no  |
        patches     | Alternate version for "boundary"      | no  |
        namedBlocks |                                       | no  |
        namedVertices |                                     | no  |
        mergeType   | Merging "points" or "topology"        | no  | topology
        checkFaceCorrespondence |                           | no  | true
        verbose     |                                       | no  | true
    \endtable

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

#include "Enum.H"
#include "block.H"
#include "PtrList.H"
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
    public PtrList<block>
{
public:

    // Typedefs

        //- The list of blocks is stored as a PtrList
        typedef PtrList<block> blockList;


    // Data Types

        //- The block merging strategy
        enum mergeStrategy
        {
            DEFAULT_MERGE,  //!< Default (TOPOLOGY), not selectable
            MERGE_TOPOLOGY, //!< "topology" merge by block topology (default)
            MERGE_POINTS    //!< "points" merge by point geometry
        };


private:

    // Static Data

        //- Names corresponding to the mergeStrategy
        static const Enum<mergeStrategy> strategyNames_;


    // Private Data

        //- Reference to mesh dictionary
        const IOdictionary& meshDict_;

        //- Output verbosity
        bool verbose_;

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

    // Static Data

        //- The default verbosity (true)
        static bool verboseOutput;


    //- Runtime type information
    ClassName("blockMesh");


    // Constructors

        //- Construct from IOdictionary for given region
        //  Default is topological merging.
        explicit blockMesh
        (
            const IOdictionary& dict,
            const word& regionName = polyMesh::defaultRegion,
            mergeStrategy strategy = mergeStrategy::DEFAULT_MERGE
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
            bool valid() const noexcept;

            //- Reference to point field defining the blockMesh.
            //  These points are \b not scaled by scaleFactor
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

            //- The points for the entire mesh.
            //  These points \b are scaled by scaleFactor
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


    // Verbosity

        //- Verbose output
        bool verbose() const noexcept;

        //- Enable/disable verbose output
        //  \return old value
        bool verbose(const bool on) noexcept;


    // Mesh Generation

        //- Create polyMesh, with cell zones
        autoPtr<polyMesh> mesh(const IOobject& io) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
