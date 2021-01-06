/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
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
    Foam::syncTools

Description
    Various tools to aid synchronizing lists across coupled patches. WIP.

    Require
    - combineOperator (e.g. sumEqOp - not sumOp!) that is defined for the
      type and combineReduce(UList\<T\>, combineOperator) should be defined.
    - null value which gets overridden by any valid value.
    - transform function

SourceFiles
    syncTools.C
    syncToolsTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef syncTools_H
#define syncTools_H

#include "Pstream.H"
#include "edgeHashes.H"
#include "bitSet.H"
#include "polyMesh.H"
#include "coupledPolyPatch.H"
#include "mapDistribute.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class polyBoundaryMesh;

/*---------------------------------------------------------------------------*\
                           Class syncTools Declaration
\*---------------------------------------------------------------------------*/

class syncTools
{
    // Private Member Functions

        //- Combine val with existing value in pointValues map at given index
        template<class T, class CombineOp>
        static void combine
        (
            Map<T>& pointValues,
            const CombineOp& cop,
            const label index,
            const T& val
        );

        //- Combine val with existing value in edgeValues at edge index
        template<class T, class CombineOp>
        static void combine
        (
            EdgeMap<T>& edgeValues,
            const CombineOp& cop,
            const edge& index,
            const T& val
        );


public:

        // Basic routines with user-supplied transformation.
        // Preferably use specialisations below.

            //- Synchronize values on selected points.
            template<class T, class CombineOp, class TransformOp>
            static void syncPointMap
            (
                const polyMesh& mesh,
                Map<T>& pointValues,
                const CombineOp& cop,
                const TransformOp& top
            );

            //- Synchronize values on selected edges.
            template<class T, class CombineOp, class TransformOp>
            static void syncEdgeMap
            (
                const polyMesh& mesh,
                EdgeMap<T>& edgeValues,
                const CombineOp& cop,
                const TransformOp& top
            );

            //- Synchronize values on all mesh points.
            template<class T, class CombineOp, class TransformOp>
            static void syncPointList
            (
                const polyMesh& mesh,
                List<T>& pointValues,
                const CombineOp& cop,
                const T& nullValue,
                const TransformOp& top
            );

            //- Synchronize values on selected mesh points.
            template<class T, class CombineOp, class TransformOp>
            static void syncPointList
            (
                const polyMesh& mesh,
                const labelUList& meshPoints,
                List<T>& pointValues,
                const CombineOp& cop,
                const T& nullValue,
                const TransformOp& top
            );

            //- Synchronize values on all mesh edges.
            template<class T, class CombineOp, class TransformOp, class FlipOp>
            static void syncEdgeList
            (
                const polyMesh& mesh,
                List<T>& edgeValues,
                const CombineOp& cop,
                const T& nullValue,
                const TransformOp& top,
                const FlipOp& fop
            );

            //- Synchronize values on selected mesh edges.
            template<class T, class CombineOp, class TransformOp, class FlipOp>
            static void syncEdgeList
            (
                const polyMesh& mesh,
                const labelList& meshEdges,
                List<T>& edgeValues,
                const CombineOp& cop,
                const T& nullValue,
                const TransformOp& top,
                const FlipOp& fop
            );

            //- Synchronize values on boundary faces only.
            template<class T, class CombineOp, class TransformOp>
            static void syncBoundaryFaceList
            (
                const polyMesh& mesh,
                UList<T>& faceValues,
                const CombineOp& cop,
                const TransformOp& top,
                const bool parRun = Pstream::parRun()
            );


        // Synchronise point-wise data

            //- Synchronize values on all mesh points.
            template<class T, class CombineOp>
            static void syncPointList
            (
                const polyMesh& mesh,
                List<T>& pointValues,
                const CombineOp& cop,
                const T& nullValue
            )
            {
                syncPointList
                (
                    mesh,
                    pointValues,
                    cop,
                    nullValue,
                    mapDistribute::transform()
                );
            }

            //- Synchronize locations on all mesh points.
            template<class CombineOp>
            static void syncPointPositions
            (
                const polyMesh& mesh,
                List<point>& positions,
                const CombineOp& cop,
                const point& nullValue
            )
            {
                syncPointList
                (
                    mesh,
                    positions,
                    cop,
                    nullValue,
                    mapDistribute::transformPosition()
                );
            }

            //- Synchronize values on selected mesh points.
            template<class T, class CombineOp>
            static void syncPointList
            (
                const polyMesh& mesh,
                const labelList& meshPoints,
                List<T>& pointValues,
                const CombineOp& cop,
                const T& nullValue
            )
            {
                syncPointList
                (
                    mesh,
                    meshPoints,
                    pointValues,
                    cop,
                    nullValue,
                    mapDistribute::transform()
                );
            }

            //- Synchronize locations on selected mesh points.
            template<class CombineOp>
            static void syncPointPositions
            (
                const polyMesh& mesh,
                const labelList& meshPoints,
                List<point>& positions,
                const CombineOp& cop,
                const point& nullValue
            )
            {
                syncPointList
                (
                    mesh,
                    meshPoints,
                    positions,
                    cop,
                    nullValue,
                    mapDistribute::transformPosition()
                );
            }


        // Synchronise edge-wise data

            //- Synchronize values on all mesh edges.
            template<class T, class CombineOp>
            static void syncEdgeList
            (
                const polyMesh& mesh,
                List<T>& edgeValues,
                const CombineOp& cop,
                const T& nullValue
            )
            {
                syncEdgeList
                (
                    mesh,
                    edgeValues,
                    cop,
                    nullValue,
                    mapDistribute::transform(),
                    noOp()
                );
            }

            //- Synchronize locations on all mesh edges.
            template<class CombineOp>
            static void syncEdgePositions
            (
                const polyMesh& mesh,
                List<point>& positions,
                const CombineOp& cop,
                const point& nullValue
            )
            {
                syncEdgeList
                (
                    mesh,
                    positions,
                    cop,
                    nullValue,
                    mapDistribute::transformPosition(),
                    noOp()
                );
            }

            //- Synchronize values on selected mesh edges.
            template<class T, class CombineOp>
            static void syncEdgeList
            (
                const polyMesh& mesh,
                const labelList& meshEdges,
                List<T>& edgeValues,
                const CombineOp& cop,
                const T& nullValue
            )
            {
                syncEdgeList
                (
                    mesh,
                    meshEdges,
                    edgeValues,
                    cop,
                    nullValue,
                    mapDistribute::transform(),
                    noOp()
                );
            }

            //- Synchronize locations on selected mesh edges.
            template<class CombineOp>
            static void syncEdgePositions
            (
                const polyMesh& mesh,
                const labelList& meshEdges,
                List<point>& positions,
                const CombineOp& cop,
                const point& nullValue
            )
            {
                syncEdgeList
                (
                    mesh,
                    meshEdges,
                    positions,
                    cop,
                    nullValue,
                    mapDistribute::transformPosition(),
                    noOp()
                );
            }



        // Synchronise face-wise data

            //- Synchronize values on boundary faces only.
            template<class T, class CombineOp>
            static void syncBoundaryFaceList
            (
                const polyMesh& mesh,
                UList<T>& faceValues,
                const CombineOp& cop
            )
            {
                syncBoundaryFaceList
                (
                    mesh,
                    faceValues,
                    cop,
                    mapDistribute::transform()
                );
            }

            //- Synchronize locations on boundary faces only.
            template<class CombineOp>
            static void syncBoundaryFacePositions
            (
                const polyMesh& mesh,
                UList<point>& positions,
                const CombineOp& cop
            )
            {
                syncBoundaryFaceList
                (
                    mesh,
                    positions,
                    cop,
                    mapDistribute::transformPosition()
                );
            }

            //- Synchronize values on all mesh faces.
            template<class T, class CombineOp>
            static void syncFaceList
            (
                const polyMesh& mesh,
                UList<T>& faceValues,
                const CombineOp& cop
            )
            {
                SubList<T> bndValues
                (
                    faceValues,
                    mesh.nBoundaryFaces(),
                    mesh.nInternalFaces()
                );

                syncBoundaryFaceList
                (
                    mesh,
                    bndValues,
                    cop,
                    mapDistribute::transform()
                );
            }

            //- Synchronize locations on all mesh faces.
            template<class CombineOp>
            static void syncFacePositions
            (
                const polyMesh& mesh,
                UList<point>& positions,
                const CombineOp& cop
            )
            {
                SubList<point> bndValues
                (
                    positions,
                    mesh.nBoundaryFaces(),
                    mesh.nInternalFaces()
                );
                syncBoundaryFaceList
                (
                    mesh,
                    bndValues,
                    cop,
                    mapDistribute::transformPosition()
                );
            }

            //- Swap coupled boundary face values. Uses eqOp
            template<class T>
            static void swapBoundaryFaceList
            (
                const polyMesh& mesh,
                UList<T>& faceValues
            )
            {
                syncBoundaryFaceList
                (
                    mesh,
                    faceValues,
                    eqOp<T>(),
                    mapDistribute::transform()
                );
            }

            //- Swap coupled positions. Uses eqOp
            static void swapBoundaryFacePositions
            (
                const polyMesh& mesh,
                UList<point>& positions
            )
            {
                syncBoundaryFaceList
                (
                    mesh,
                    positions,
                    eqOp<point>(),
                    mapDistribute::transformPosition()
                );
            }

            //- Swap coupled face values. Uses eqOp
            template<class T>
            static void swapFaceList
            (
                const polyMesh& mesh,
                UList<T>& faceValues
            )
            {
                SubList<T> bndValues
                (
                    faceValues,
                    mesh.nBoundaryFaces(),
                    mesh.nInternalFaces()
                );
                syncBoundaryFaceList
                (
                    mesh,
                    bndValues,
                    eqOp<T>(),
                    mapDistribute::transform()
                );
            }

            //- Swap to obtain neighbour cell values for all boundary faces
            template<class T>
            static void swapBoundaryCellList
            (
                const polyMesh& mesh,
                const UList<T>& cellData,
                List<T>& neighbourCellData
            );

            //- Swap to obtain neighbour cell positions for all boundary faces
            static void swapBoundaryCellPositions
            (
                const polyMesh& mesh,
                const UList<point>& cellData,
                List<point>& neighbourCellData
            );


        // Sparse versions

            //- Synchronize values on selected points.
            template<class T, class CombineOp>
            static void syncPointMap
            (
                const polyMesh& mesh,
                Map<T>& pointValues,
                const CombineOp& cop
            )
            {
                syncPointMap
                (
                    mesh,
                    pointValues,
                    cop,
                    mapDistribute::transform()
                );
            }

            //- Synchronize locations on selected points.
            template<class CombineOp>
            static void syncPointPositions
            (
                const polyMesh& mesh,
                Map<point>& positions,
                const CombineOp& cop
            )
            {
                syncPointMap
                (
                    mesh,
                    positions,
                    cop,
                    mapDistribute::transformPosition()
                );
            }

            //- Synchronize values on selected edges. Edges are represented
            //  by the two vertices that make it up so global edges never get
            //  constructed.
            template<class T, class CombineOp>
            static void syncEdgeMap
            (
                const polyMesh& mesh,
                EdgeMap<T>& edgeValues,
                const CombineOp& cop
            )
            {
                syncEdgeMap
                (
                    mesh,
                    edgeValues,
                    cop,
                    mapDistribute::transform()
                );
            }

            //- Synchronize locations on selected edges.
            template<class CombineOp>
            static void syncEdgePositions
            (
                const polyMesh& mesh,
                EdgeMap<point>& edgePositions,
                const CombineOp& cop
            )
            {
                syncEdgeMap
                (
                    mesh,
                    edgePositions,
                    cop,
                    mapDistribute::transformPosition()
                );
            }


        // PackedList versions

            //- Synchronize face values from PackedList/bitSet
            //
            //  \param mesh The mesh
            //  \param isBoundaryOnly True if faceValues are for the boundary
            //      only and not the entire mesh. This determines the face
            //      offset when accessing values.
            //  \param faceValues The face values to synchronize
            //  \param cop The combine operation
            //  \param parRun True if this is a parallel simulation
            template<unsigned Width, class CombineOp>
            static void syncFaceList
            (
                const polyMesh& mesh,
                const bool isBoundaryOnly,
                PackedList<Width>& faceValues,
                const CombineOp& cop,
                const bool parRun = Pstream::parRun()
            );

            //- Synchronize mesh face values from PackedList/bitSet
            template<unsigned Width, class CombineOp>
            static void syncFaceList
            (
                const polyMesh& mesh,
                PackedList<Width>& faceValues,
                const CombineOp& cop,
                const bool parRun = Pstream::parRun()
            );

            //- Synchronize boundary face values from PackedList/bitSet
            template<unsigned Width, class CombineOp>
            static void syncBoundaryFaceList
            (
                const polyMesh& mesh,
                PackedList<Width>& faceValues,
                const CombineOp& cop,
                const bool parRun = Pstream::parRun()
            );

            //- Swap coupled face values. Uses eqOp
            template<unsigned Width>
            static void swapFaceList
            (
                const polyMesh& mesh,
                PackedList<Width>& faceValues
            );

            //- Swap coupled boundary face values. Uses eqOp
            template<unsigned Width>
            static void swapBoundaryFaceList
            (
                const polyMesh& mesh,
                PackedList<Width>& faceValues
            );

            template<unsigned Width, class CombineOp>
            static void syncPointList
            (
                const polyMesh& mesh,
                PackedList<Width>& pointValues,
                const CombineOp& cop,
                const unsigned int nullValue
            );

            template<unsigned Width, class CombineOp>
            static void syncEdgeList
            (
                const polyMesh& mesh,
                PackedList<Width>& edgeValues,
                const CombineOp& cop,
                const unsigned int nullValue
            );


        // Other

            //- Get per point whether it is uncoupled or a master of a
            //- coupled set of points
            static bitSet getMasterPoints(const polyMesh& mesh);

            //- Get per edge whether it is uncoupled or a master of a
            //- coupled set of edges
            static bitSet getMasterEdges(const polyMesh& mesh);

            //- Get per face whether it is uncoupled or a master of a
            //- coupled set of faces
            static bitSet getMasterFaces(const polyMesh& mesh);

            //- Get per face whether it is internal or a master of a
            //- coupled set of faces
            static bitSet getInternalOrMasterFaces(const polyMesh& mesh);

            //- Get per face whether it is internal or coupled
            static bitSet getInternalOrCoupledFaces(const polyMesh& mesh);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "syncToolsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
