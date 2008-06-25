/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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
    Foam::autoHexMeshDriver

Description
    main meshing driver.

SourceFiles
    autoHexMeshDriver.C
    autoHexMeshDriverSnap.C
    autoHexMeshDriverLayers.C
    autoHexMeshDriverShrink.C
    autoHexMeshDriverTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef autoHexMeshDriver_H
#define autoHexMeshDriver_H

#include "autoPtr.H"
#include "dictionary.H"
#include "pointField.H"
#include "boolList.H"
#include "Switch.H"
#include "PackedList.H"
#include "wallPoint.H"
#include "indirectPrimitivePatch.H"
#include "featureEdgeMesh.H"
#include "searchableSurface.H"
#include "refinementSurfaces.H"
#include "meshRefinement.H"
#include "decompositionMethod.H"
#include "fvMeshDistribute.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class fvMesh;
class pointMesh;
class motionSmoother;
class removePoints;
class pointSet;
class pointData;
class faceSet;
class addPatchCellLayer;
class mapDistributePolyMesh;

/*---------------------------------------------------------------------------*\
                           Class autoHexMeshDriver Declaration
\*---------------------------------------------------------------------------*/

class autoHexMeshDriver
{
    // Static data members

        //- Default angle for faces to be convcave
        static const scalar defaultConcaveAngle;

        //- Extrusion controls
        enum extrudeMode
        {
            NOEXTRUDE,      /*!< Do not extrude. No layers added. */
            EXTRUDE,        /*!< Extrude */
            EXTRUDEREMOVE   /*!< Extrude but afterwards remove added faces locally */
        };


    // Private classes

        //- Combine operator class for equalizing displacements.
        class minMagEqOp
        {
        public:

            void operator()(vector& x, const vector& y) const
            {
                if (magSqr(y) < magSqr(x))
                {
                    x = y;
                }
            }
        };

        //- Combine operator class to combine normal with other normal.
        class nomalsCombine
        {
        public:

            void operator()(vector& x, const vector& y) const
            {
                if (y != wallPoint::greatPoint)
                {
                    if (x == wallPoint::greatPoint)
                    {
                        x = y;
                    }
                    else
                    {
                        x *= (x&y);
                    }
                }
            }
        };



    // Private data

        //- Reference to mesh
        fvMesh& mesh_;

        //- Input dictionary
        const dictionary dict_;

        //- Debug level
        const label debug_;

        //- Total number of cells
        const label maxGlobalCells_;

        //- Per processor max number of cells
        const label maxLocalCells_;

        //- When to stop refining
        const label minRefineCells_;

        //- Curvature
        const scalar curvature_;

        //- Number of layers between different refinement levels
        const label nBufferLayers_;

        //- Areas to keep
        const pointField keepPoints_;

        //- Merge distance
        const scalar mergeDist_;



        //- Explicit features
        PtrList<featureEdgeMesh> featureMeshes_;
        //- Per feature the refinement level
        labelList featureLevels_;

        //- Shells
        PtrList<searchableSurface> shells_;
        //- Per shell the refinement level
        labelList shellLevels_;
        //- Per shell whether to refine inside or outside
        boolList shellRefineInside_;

        //- Surfaces with refinement levels built-in
        autoPtr<refinementSurfaces> surfacesPtr_;

        //- Per refinement surface region the patch
        labelList globalToPatch_;

        //- Mesh refinement engine
        autoPtr<meshRefinement> meshRefinerPtr_;

        //- Decomposition engine
        autoPtr<decompositionMethod> decomposerPtr_;

        //- Mesh distribution engine
        autoPtr<fvMeshDistribute> distributorPtr_;



    // Private Member Functions

        // Refinement

            //- Calculate merge distance. Check against writing tolerance.
            scalar getMergeDistance(const scalar mergeTol) const;

            // Return per keeppoint -1 or the local cell label the point is in.
            // Guaranteed to be only on one processor.
            labelList findCells(const pointField&) const;

            static void orientOutside(PtrList<searchableSurface>&);

            //- Read feature edges
            label readFeatureEdges(const PtrList<dictionary>&);

        // Snapping

             //- Split surfaces into non-zoned and zones ones
            void getZonedSurfaces(labelList&, labelList&) const;

            //- Get faces to repatch. Returns map from face to patch.
            Map<label> getZoneBafflePatches(const bool allowBoundary) const;

            //- Calculates (geometric) shared points
            static label getCollocatedPoints
            (
                const scalar tol,
                const pointField&,
                PackedList<1>&
            );

            //- Calculate displacement per patch point to smooth out patch.
            //  Quite complicated in determining which points to move where.
            pointField smoothPatchDisplacement(const motionSmoother&) const;

            //- Check that face zones are synced
            void checkCoupledFaceZones() const;

            //- Per edge distance to patch
            static tmp<scalarField> edgePatchDist
            (
                const pointMesh&,
                const indirectPrimitivePatch&
            );

            //- Write displacement as .obj file.
            static void dumpMove
            (
                const fileName&,
                const pointField&,
                const pointField&
            );

            //- Check displacement is outwards pointing
            static bool outwardsDisplacement
            (
                const indirectPrimitivePatch&,
                const vectorField&
            );

        // Face merging

            //- Merge patch faces. Undo until no checkMesh errors.
            label mergePatchFacesUndo
            (
                const scalar minCos,
                const scalar concaveCos,
                const dictionary&
            );

            //- Remove points.
            autoPtr<mapPolyMesh> doRemovePoints
            (
                removePoints& pointRemover,
                const boolList& pointCanBeDeleted
            );

            //- Restore faces (which contain removed points)
            autoPtr<mapPolyMesh> doRestorePoints
            (
                removePoints& pointRemover,
                const labelList& facesToRestore
            );

            //- Return candidateFaces that are also in set.
            labelList collectFaces
            (
                const labelList& candidateFaces,
                const labelHashSet& set
            ) const;

            //- Pick up faces of cells of faces in set.
            labelList growFaceCellFace(const labelHashSet&) const;

            //- Remove points not used by any face or points used by only
            //  two faces where the edges are in line
            label mergeEdgesUndo(const scalar minCos, const dictionary&);


        // Layers

            //- For debugging: Dump displacement to .obj files
            static void dumpDisplacement
            (
                const fileName&,
                const indirectPrimitivePatch&,
                const vectorField&,
                const List<extrudeMode>&
            );

            //- Check that primitivePatch is not multiply connected.
            //  Collect non-manifold points in pointSet.
            static void checkManifold
            (
                const indirectPrimitivePatch&,
                pointSet& nonManifoldPoints
            );


            // Static extrusion setup

                //- Unset extrusion on point. Returns true if anything unset.
                static bool unmarkExtrusion
                (
                    const label patchPointI,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                );

                //- Unset extrusion on face. Returns true if anything unset.
                static bool unmarkExtrusion
                (
                    const face& localFace,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                );

                //- No extrusion at non-manifold points.
                void handleNonManifolds
                (
                    const indirectPrimitivePatch& pp,
                    const labelList& meshEdges,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                ) const;

                //- No extrusion on feature edges. Assumes non-manifold
                //  edges already handled.
                void handleFeatureAngle
                (
                    const indirectPrimitivePatch& pp,
                    const labelList& meshEdges,
                    const scalar minCos,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                ) const;

                //- No extrusion on warped faces
                void handleWarpedFaces
                (
                    const indirectPrimitivePatch& pp,
                    const scalar faceRatio,
                    const scalar edge0Len,
                    const labelList& cellLevel,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                ) const;

                //- Determine the number of layers per point from the number of
                //  layers per surface.
                void setNumLayers
                (
                    const labelList& patchIDs,
                    const indirectPrimitivePatch& pp,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                ) const;

                //- Grow no-extrusion layer.
                static void growNoExtrusion
                (
                    const indirectPrimitivePatch& pp,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                );

                //- Calculate pointwise wanted and minimum thickness.
                //  thickness: wanted thickness
                //  minthickness: when to give up and not extrude
                void calculateLayerThickness
                (
                    const scalar expansionRatio,
                    const scalar finalLayerRatio,
                    const scalar relMinThickness,
                    const indirectPrimitivePatch& pp,
                    const labelList& cellLevel,
                    const labelList& patchNLayers,
                    const scalar edge0Len,
                    scalarField& thickness,
                    scalarField& minThickness
                ) const;


            // Extrusion execution

                //- Synchronize displacement among coupled patches.
                void syncPatchDisplacement
                (
                    const motionSmoother& meshMover,
                    const scalarField& minThickness,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                ) const;

                //- Get nearest point on surface to snap to
                void getPatchDisplacement
                (
                    const motionSmoother& meshMover,
                    const scalarField& thickness,
                    const scalarField& minThickness,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                ) const;

                //- Truncates displacement
                // - for all patchFaces in the faceset displacement gets set
                //   to zero
                // - all displacement < minThickness gets set to zero
                label truncateDisplacement
                (
                    const motionSmoother& meshMover,
                    const scalarField& minThickness,
                    const faceSet& illegalPatchFaces,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                ) const;

                //- Setup layer information (at points and faces) to
                //  modify mesh topology in
                //  regions where layer mesh terminates. Guarantees an
                //  optional slow decreasing of the number of layers.
                //  Returns the number of layers per face and per point
                //  to go into the actual layer addition engine.
                void setupLayerInfoTruncation
                (
                    const motionSmoother& meshMover,
                    const labelList& patchNLayers,
                    const List<extrudeMode>& extrudeStatus,
                    const label nBufferCellsNoExtrude,
                    labelList& nPatchPointLayers,
                    labelList& nPatchFaceLayers
                ) const;

                //- Does any of the cells use a face from faces?
                static bool cellsUseFace
                (
                    const polyMesh& mesh,
                    const labelList& cellLabels,
                    const labelHashSet& faces
                );

                //- Checks the newly added cells and locally unmarks points
                //  so they will not get extruded next time round. Returns
                //  global number of unmarked points (0 if all was fine)
                static label checkAndUnmark
                (
                    const addPatchCellLayer& addLayer,
                    const dictionary& motionDict,
                    const indirectPrimitivePatch& pp,
                    const fvMesh&,

                    pointField& patchDisp,
                    labelList& patchNLayers,
                    List<extrudeMode>& extrudeStatus
                );

                //- Count global number of extruded faces
                static label countExtrusion
                (
                    const indirectPrimitivePatch& pp,
                    const List<extrudeMode>& extrudeStatus
                );

                //- Collect layer faces and layer cells into bools
                //  for ease of handling
                static void getLayerCellsFaces
                (
                    const polyMesh&,
                    const addPatchCellLayer&,
                    boolList&,
                    boolList&
                );

            // Mesh shrinking (to create space for layers)

                //- Average field (over all subset of mesh points) by
                //  summing contribution from edges. Global parallel since only
                //  does master edges for coupled edges.
                template<class Type>
                void averageNeighbours
                (
                    const PackedList<1>& isMasterEdge,
                    const labelList& meshEdges,
                    const labelList& meshPoints,
                    const edgeList& edges,
                    const scalarField& invSumWeight,
                    const Field<Type>& data,
                    Field<Type>& average
                ) const;

                //- Calculate inverse sum of edge weights (currently always 1.0)
                void sumWeights
                (
                    const PackedList<1>& isMasterEdge,
                    const labelList& meshEdges,
                    const labelList& meshPoints,
                    const edgeList& edges,
                    scalarField& invSumWeight
                ) const;

                //- Smooth scalar field on patch
                void smoothField
                (
                    const motionSmoother& meshMover,
                    const PackedList<1>& isMasterEdge,
                    const labelList& meshEdges,
                    const scalarField& fieldMin,
                    const label& nSmoothDisp,
                    scalarField& field
                ) const;

                //- Smooth normals on patch.
                void smoothPatchNormals
                (
                    const motionSmoother& meshMover,
                    const PackedList<1>& isMasterEdge,
                    const labelList& meshEdges,
                    const label nSmoothDisp,
                    pointField& normals
                ) const;

                //- Smooth normals in interior.
                void smoothNormals
                (
                    const label nSmoothDisp,
                    const PackedList<1>& isMasterEdge,
                    const labelList& fixedPoints,
                    pointVectorField& normals
                ) const;

                bool isMaxEdge
                (
                    const List<pointData>&,
                    const label edgeI,
                    const scalar minCos
                ) const;

                //- Stop layer growth where mesh wraps around edge with a
                //  large feature angle
                void handleFeatureAngleLayerTerminations
                (
                    const indirectPrimitivePatch& pp,
                    const scalar minCos,
                    List<extrudeMode>& extrudeStatus,
                    pointField& patchDisp,
                    labelList& patchNLayers,
                    label& nPointCounter
                ) const;

                //- Find isolated islands (points, edges and faces and
                // layer terminations)
                // in the layer mesh and stop any layer growth at these points.
                void findIsolatedRegions
                (
                    const indirectPrimitivePatch& pp,
                    const PackedList<1>& isMasterEdge,
                    const labelList& meshEdges,
                    const scalar minCosLayerTermination,
                    scalarField& field,
                    List<extrudeMode>& extrudeStatus,
                    pointField& patchDisp,
                    labelList& patchNLayers
                ) const;

                // Calculate medial axis fields
                void medialAxisSmoothingInfo
                (
                    const motionSmoother& meshMover,
                    const label nSmoothNormals,
                    const label nSmoothSurfaceNormals,
                    const scalar minMedianAxisAngleCos,

                    pointVectorField& dispVec,
                    pointScalarField& medialRatio,
                    pointScalarField& medialDist
                ) const;

                //- Main routine to shrink mesh
                void shrinkMeshMedialDistance
                (
                    motionSmoother& meshMover,
                    const label nSmoothThickness,
                    const scalar maxThicknessToMedialRatio,
                    const label nAllowableErrors,
                    const label nSnap,
                    const scalar minCosLayerTermination,

                    const scalarField& layerThickness,
                    const scalarField& minThickness,

                    const pointVectorField& dispVec,
                    const pointScalarField& medialRatio,
                    const pointScalarField& medialDist,

                    List<extrudeMode>& extrudeStatus,
                    pointField& patchDisp,
                    labelList& patchNLayers
                ) const;


        //- Disallow default bitwise copy construct
        autoHexMeshDriver(const autoHexMeshDriver&);

        //- Disallow default bitwise assignment
        void operator=(const autoHexMeshDriver&);

public:

    //- Runtime type information
    ClassName("autoHexMeshDriver");


    // Constructors

        //- Construct from dictionary and mesh to modify
        autoHexMeshDriver
        (
            fvMesh& mesh,
            const dictionary& meshDict,
            const dictionary& decomposeDict
        );

        //- Alternative constructor from top-level controldictionary and
        //  refinement specific dictionary
        autoHexMeshDriver
        (
            fvMesh& mesh,
            const dictionary& controlDict,
            const dictionary& refineDict,
            const dictionary& decomposeDict
        );


    // Member Functions

        // Access

            //- reference to mesh
            const fvMesh& mesh() const
            {
                return mesh_;
            }
            fvMesh& mesh()
            {
                return mesh_;
            }

            //- Surfaces to base refinement on
            const refinementSurfaces& surfaces() const
            {
                return surfacesPtr_();
            }

            //- Per refinementsurface, per region the patch
            const labelList& globalToPatch() const
            {
                return globalToPatch_;
            }


        // Refinement

            //- Refine around explicit feature edges
            label featureEdgeRefine
            (
                const PtrList<dictionary>& featDicts,
                const label maxIter,
                const label minRefine
            );

            //- Refine at surface intersections
            label surfaceOnlyRefine(const label maxIter);

            //- Remove cells not reachable from keep points
            void removeInsideCells(const label nBufferLayers);

            //- Refine volume regions (interior of shells)
            label shellRefine(const label maxIter);

            //- Introduce baffles; remove unreachable mesh parts
            //  handleSnapProblems : whether to remove free floating cells
            void baffleAndSplitMesh(const bool handleSnapProblems);

            //- Move cells to zones
            void zonify();

            //- Split and recombine baffles to get rid of single face baffles.
            void splitAndMergeBaffles(const bool handleSnapProblems);

            //- Merge multiple boundary faces on single cell
            void mergePatchFaces();

            //- Redecompose according to cell count
            //  keepZoneFaces : find all faceZones from zoned surfaces and keep
            //                  owner and neighbour together
            //  keepBaffles   : find all baffles and keep them together
            autoPtr<mapDistributePolyMesh> balance
            (
                const bool keepZoneFaces,
                const bool keepBaffles
            );

            //- Write mesh
            void writeMesh(const string&) const;


        // Snapping

            //- Create baffles for faces straddling zoned surfaces. Return
            //  baffles.
            autoPtr<mapPolyMesh> createZoneBaffles(List<labelPair>&);

            //- Merge baffles.
            autoPtr<mapPolyMesh> mergeZoneBaffles(const List<labelPair>&);

            //- Calculate edge length per patch point.
            scalarField calcSnapDistance
            (
                const dictionary& snapDict,
                const indirectPrimitivePatch&
            ) const;

            //- Get patches generated for surfaces.
            labelList getSurfacePatches() const;

            //- Smooth the mesh (patch and internal) to increase visibility
            //  of surface points (on castellated mesh) w.r.t. surface.
            void preSmoothPatch
            (
                const dictionary& snapDict,
                const label nInitErrors,
                const List<labelPair>& baffles,
                motionSmoother&
            ) const;

            //- Per patch point calculate point on nearest surface. Set as
            //  boundary conditions of motionSmoother displacement field. Return
            //  displacement of patch points.
            vectorField calcNearestSurface
            (
                const scalarField& snapDist,
                motionSmoother& meshMover
            ) const;

            //- Smooth the displacement field to the internal.
            void smoothDisplacement
            (
                const dictionary& snapDict,
                motionSmoother&
            ) const;

            //- Do the hard work: move the mesh according to displacement,
            //  locally relax the displacement.
            void scaleMesh
            (
                const dictionary& snapDict,
                const label nInitErrors,
                const List<labelPair>& baffles,
                motionSmoother&
            );


        // Layers

            //- Merge patch faces on same cell.
            void mergePatchFacesUndo
            (
                const dictionary& shrinkDict,
                const dictionary& motionDict
            );

            //- Check that mesh outside is not multiply connected.
            void checkMeshManifold() const;

            //- Add cell layers
            void addLayers
            (
                const dictionary& shrinkDict,
                const dictionary& motionDict,
                const label nAllowableErrors,
                motionSmoother&
            );


        // Other

            //- Do all refinement.
            void doRefine
            (
                const dictionary& refineDict,
                const bool prepareForSnapping
            );

            //- Do all snapping.
            void doSnap
            (
                const dictionary& snapDict,
                const dictionary& motionDict
            );

            //- Do alllayer addition.
            void doLayers
            (
                const dictionary& shrinkDict,
                const dictionary& motionDict
            );

            //- Do all : refine, snap, layers
            void doMesh();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "autoHexMeshDriverTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
