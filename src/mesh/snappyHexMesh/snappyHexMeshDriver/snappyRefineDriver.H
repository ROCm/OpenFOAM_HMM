/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd.
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
    Foam::snappyRefineDriver

Description

SourceFiles
    snappyRefineDriver.C

\*---------------------------------------------------------------------------*/

#ifndef snappyRefineDriver_H
#define snappyRefineDriver_H

#include "wordPairHashTable.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class refinementParameters;
class snapParameters;

class meshRefinement;
class decompositionMethod;
class fvMeshDistribute;

/*---------------------------------------------------------------------------*\
                           Class snappyRefineDriver Declaration
\*---------------------------------------------------------------------------*/

class snappyRefineDriver
{
    // Private data

        //- Mesh+surface
        meshRefinement& meshRefiner_;

        //- Reference to decomposition method
        decompositionMethod& decomposer_;

        //- Reference to mesh distribution engine
        fvMeshDistribute& distributor_;

        //- From surface region to patch
        const labelList globalToMasterPatch_;

        //- From surface region to patch
        const labelList globalToSlavePatch_;


    // Private Member Functions

        //- Refine all cells pierced by explicit feature edges
        label featureEdgeRefine
        (
            const refinementParameters& refineParams,
            const label maxIter,
            const label minRefine
        );

        //- Refine all cells containing small surface features
        label smallFeatureRefine
        (
            const refinementParameters& refineParams,
            const label maxIter
        );

        //- Refine all cells interacting with the surface
        label surfaceOnlyRefine
        (
            const refinementParameters& refineParams,
            const label maxIter
        );

        //- Refine all cells in small gaps
        label gapOnlyRefine
        (
            const refinementParameters& refineParams,
            const label maxIter
        );

        //- Refine all cells in large gaps
        label bigGapOnlyRefine
        (
            const refinementParameters& refineParams,
            const bool spreadGapSize,
            const label maxIter
        );

        //- Refine cells with almost all sides refined
        label danglingCellRefine
        (
            const refinementParameters& refineParams,
            const label nFaces,
            const label maxIter
        );

        //- Refine cells with opposite faces with differing refinement level
        label refinementInterfaceRefine
        (
            const refinementParameters& refineParams,
            const label maxIter
        );

        //- Remove all cells within intersected region
        void removeInsideCells
        (
            const refinementParameters& refineParams,
            const label nBufferLayers
        );

        //- Remove all cells inside/outside shell
        label shellRefine
        (
            const refinementParameters& refineParams,
            const label maxIter
        );

        //- Add baffles and remove unreachable cells
        void baffleAndSplitMesh
        (
            const refinementParameters& refineParams,
            const snapParameters& snapParams,
            const bool handleSnapProblems,
            const dictionary& motionDict
        );

        //- Add zones
        void zonify
        (
            const refinementParameters& refineParams,
            wordPairHashTable& zonesToFaceZone
        );

        void splitAndMergeBaffles
        (
            const refinementParameters& refineParams,
            const snapParameters& snapParams,
            const bool handleSnapProblems,
            const dictionary& motionDict
        );

        //- Merge refined boundary faces (from exposing coarser cell)
        void mergePatchFaces
        (
            const bool geometricMerge,
            const refinementParameters& refineParams,
            const dictionary& motionDict
        );

        //- Disallow default bitwise copy construct
        snappyRefineDriver(const snappyRefineDriver&);

        //- Disallow default bitwise assignment
        void operator=(const snappyRefineDriver&);


public:

    //- Runtime type information
    ClassName("snappyRefineDriver");


    // Constructors

        //- Construct from components
        snappyRefineDriver
        (
            meshRefinement& meshRefiner,
            decompositionMethod& decomposer,
            fvMeshDistribute& distributor,
            const labelList& globalToMasterPatch,
            const labelList& globalToSlavePatch
        );


    // Member Functions

        //- Do all the refinement
        void doRefine
        (
            const dictionary& refineDict,
            const refinementParameters& refineParams,
            const snapParameters& snapParams,
            const bool prepareForSnapping,
            const bool mergePatchFaces,
            const dictionary& motionDict
        );

        //- Helper: add faceZones and patches
        static void addFaceZones
        (
            meshRefinement& meshRefiner,
            const refinementParameters& refineParams,
            const HashTable<Pair<word>>& faceZoneToPatches
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
