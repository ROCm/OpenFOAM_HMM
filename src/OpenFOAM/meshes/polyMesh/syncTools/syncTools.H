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
    Foam::syncTools

Description
    Various tools to aid synchronizing lists across coupled patches.

    Require
    - combineOperator (e.g. sumEqOp - not sumOp!) that is defined for the
      type and combineReduce(UList\<T\>, combineOperator) should be defined.
    - null value which gets overridden by any valid value.
    - transform function

    Can apply coordinate rotation/separation on cyclics but only for faces
    or if there is a single rotation/separation tensor.

SourceFiles
    syncTools.C
    syncToolsTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef syncTools_H
#define syncTools_H

#include "UList.H"
#include "Pstream.H"
#include "Map.H"
#include "EdgeMap.H"
#include "PackedBoolList.H"
#include "polyMesh.H"
#include "coupledPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyBoundaryMesh;

/*---------------------------------------------------------------------------*\
                           Class syncTools Declaration
\*---------------------------------------------------------------------------*/

class syncTools
{
    // Private Member Functions

        //- Check whether uses couples.
        static bool hasCouples(const polyBoundaryMesh&);

        //- Combine value with existing value in map.
        template <class T, class CombineOp>
        static void combine
        (
            Map<T>& pointValues,
            const CombineOp& cop,
            const label index,
            const T& val
        );

        //- Combine val with existing value at (implicit index) e.
        template <class T, class CombineOp>
        static void combine
        (
            EdgeMap<T>& edgeValues,
            const CombineOp& cop,
            const edge& index,
            const T& val
        );


public:

    // Public classes

        class transform
        {
        public:
            //- Transform Field
            template<class T>
            void operator()(const coupledPolyPatch& cpp, Field<T>& fld) const
            {
                cpp.transform(fld);
            }

            //- Transform generic container
            template<class T, template<class> class Container>
            void operator()(const coupledPolyPatch& cpp, Container<T>& map)
            const
            {
                Field<T> fld(map.size());
                label i = 0;
                forAllConstIter(typename Container<T>, map, iter)
                {
                    fld[i++] = iter();
                }
                cpp.transform(fld);
                i = 0;
                forAllIter(typename Container<T>, map, iter)
                {
                    iter() = fld[i++];
                }
            }
        };

        class transformPosition
        {
        public:
            void operator()(const coupledPolyPatch& cpp, pointField& fld) const
            {
                cpp.transformPosition(fld);
            }
            template<template<class> class Container>
            void operator()(const coupledPolyPatch& cpp, Container<point>& map)
            const
            {
                Field<point> fld(map.size());
                label i = 0;
                forAllConstIter(typename Container<point>, map, iter)
                {
                    fld[i++] = iter();
                }
                cpp.transform(fld);
                i = 0;
                forAllIter(typename Container<point>, map, iter)
                {
                    iter() = fld[i++];
                }
            }
        };


        // Basic routines with user-supplied transformation. Preferably
        // use specialisations below.

            //- Synchronize values on selected points.
            template <class T, class CombineOp, class TransformOp>
            static void syncPointMap
            (
                const polyMesh&,
                Map<T>& pointValues,
                const CombineOp& cop,
                const TransformOp& top
            );

            //- Synchronize values on selected edges.
            template <class T, class CombineOp, class TransformOp>
            static void syncEdgeMap
            (
                const polyMesh&,
                EdgeMap<T>& edgeValues,
                const CombineOp& cop,
                const TransformOp& top
            );

            //- Synchronize values on all mesh points.
            template <class T, class CombineOp, class TransformOp>
            static void syncPointList
            (
                const polyMesh&,
                UList<T>&,
                const CombineOp& cop,
                const T& nullValue,
                const TransformOp& top
            );

            //- Synchronize values on selected mesh points.
            template <class T, class CombineOp, class TransformOp>
            static void syncPointList
            (
                const polyMesh&,
                const labelList& meshPoints,
                UList<T>&,
                const CombineOp& cop,
                const T& nullValue,
                const TransformOp& top
            );

            //- Synchronize values on all mesh edges.
            template <class T, class CombineOp, class TransformOp>
            static void syncEdgeList
            (
                const polyMesh&,
                UList<T>&,
                const CombineOp& cop,
                const T& nullValue,
                const TransformOp& top
            );

            //- Synchronize values on boundary faces only.
            template <class T, class CombineOp, class TransformOp>
            static void syncBoundaryFaceList
            (
                const polyMesh&,
                UList<T>&,
                const CombineOp& cop,
                const TransformOp& top
            );


        // Synchronise point-wise data

            //- Synchronize values on all mesh points.
            template <class T, class CombineOp>
            static void syncPointList
            (
                const polyMesh& mesh,
                UList<T>& l,
                const CombineOp& cop,
                const T& nullValue
            )
            {
                syncPointList(mesh, l, cop, nullValue, transform());
            }

            //- Synchronize locations on all mesh points.
            template <class CombineOp>
            static void syncPointPositions
            (
                const polyMesh& mesh,
                UList<point>& l,
                const CombineOp& cop,
                const point& nullValue
            )
            {
                syncPointList(mesh, l, cop, nullValue, transformPosition());
            }

            //- Synchronize values on selected mesh points.
            template <class T, class CombineOp>
            static void syncPointList
            (
                const polyMesh& mesh,
                const labelList& meshPoints,
                UList<T>& l,
                const CombineOp& cop,
                const T& nullValue
            )
            {
                syncPointList
                (
                    mesh,
                    meshPoints,
                    l,
                    cop,
                    nullValue,
                    transform()
                );
            }

            //- Synchronize locations on selected mesh points.
            template <class CombineOp>
            static void syncPointPositions
            (
                const polyMesh& mesh,
                const labelList& meshPoints,
                UList<point>& l,
                const CombineOp& cop,
                const point& nullValue
            )
            {
                syncPointList
                (
                    mesh,
                    meshPoints,
                    l,
                    cop,
                    nullValue,
                    transformPosition()
                );
            }


        // Synchronise edge-wise data

            //- Synchronize values on all mesh edges.
            template <class T, class CombineOp>
            static void syncEdgeList
            (
                const polyMesh& mesh,
                UList<T>& l,
                const CombineOp& cop,
                const T& nullValue
            )
            {
                syncEdgeList(mesh, l, cop, nullValue, transform());
            }

            //- Synchronize values on all mesh edges.
            template <class CombineOp>
            static void syncEdgePositions
            (
                const polyMesh& mesh,
                UList<point>& l,
                const CombineOp& cop,
                const point& nullValue
            )
            {
                syncEdgeList(mesh, l, cop, nullValue, transformPosition());
            }


        // Synchronise face-wise data

            //- Synchronize values on boundary faces only.
            template <class T, class CombineOp>
            static void syncBoundaryFaceList
            (
                const polyMesh& mesh,
                UList<T>& l,
                const CombineOp& cop
            )
            {
                syncBoundaryFaceList(mesh, l, cop, transform());
            }

            //- Synchronize locations on boundary faces only.
            template <class CombineOp>
            static void syncBoundaryFacePositions
            (
                const polyMesh& mesh,
                UList<point>& l,
                const CombineOp& cop
            )
            {
                syncBoundaryFaceList(mesh, l, cop, transformPosition());
            }

            //- Synchronize values on all mesh faces.
            template <class T, class CombineOp>
            static void syncFaceList
            (
                const polyMesh& mesh,
                UList<T>& l,
                const CombineOp& cop
            )
            {
                SubList<T> bndValues
                (
                    l,
                    mesh.nFaces()-mesh.nInternalFaces(),
                    mesh.nInternalFaces()
                );

                syncBoundaryFaceList(mesh, bndValues, cop, transform());
            }

            //- Synchronize locations on all mesh faces.
            template <class CombineOp>
            static void syncFacePositions
            (
                const polyMesh& mesh,
                UList<point>& l,
                const CombineOp& cop
            )
            {
                SubList<point> bndValues
                (
                    l,
                    mesh.nFaces()-mesh.nInternalFaces(),
                    mesh.nInternalFaces()
                );
                syncBoundaryFaceList(mesh, bndValues, cop, transformPosition());
            }

            //- Swap coupled boundary face values.
            template <class T>
            static void swapBoundaryFaceList
            (
                const polyMesh& mesh,
                UList<T>& l
            )
            {
                syncBoundaryFaceList(mesh, l, eqOp<T>(), transform());
            }

             //- Swap coupled positions.
            template <class T>
            static void swapBoundaryFacePositions
            (
                const polyMesh& mesh,
                UList<T>& l
            )
            {
                syncBoundaryFaceList(mesh, l, eqOp<T>(), transformPosition());
            }

            //- Swap coupled face values.
            template <class T>
            static void swapFaceList
            (
                const polyMesh& mesh,
                UList<T>& l
            )
            {
                SubList<T> bndValues
                (
                    l,
                    mesh.nFaces()-mesh.nInternalFaces(),
                    mesh.nInternalFaces()
                );
                syncBoundaryFaceList(mesh, bndValues, eqOp<T>(), transform());
            }

        // Sparse versions

            //- Synchronize values on selected points.
            template <class T, class CombineOp>
            static void syncPointMap
            (
                const polyMesh& mesh,
                Map<T>& l,
                const CombineOp& cop
            )
            {
                syncPointMap(mesh, l, cop, transform());
            }

            //- Synchronize locations on selected points.
            template <class CombineOp>
            static void syncPointPositions
            (
                const polyMesh& mesh,
                Map<point>& l,
                const CombineOp& cop
            )
            {
                syncPointMap(mesh, l, cop, transformPosition());
            }

            //- Synchronize values on selected edges. Edges are represented
            //  by the two vertices that make it up so global edges never get
            //  constructed.
            template <class T, class CombineOp>
            static void syncEdgeMap
            (
                const polyMesh& mesh,
                EdgeMap<T>& l,
                const CombineOp& cop
            )
            {
                syncEdgeMap(mesh, l, cop, transform());
            }

            //- Synchronize locations on selected edges.
            template <class T, class CombineOp>
            static void syncEdgePositions
            (
                const polyMesh& mesh,
                EdgeMap<T>& l,
                const CombineOp& cop
            )
            {
                syncEdgeMap(mesh, l, cop, transformPosition());
            }

        // PackedList versions

            template <unsigned nBits, class CombineOp>
            static void syncFaceList
            (
                const polyMesh& mesh,
                PackedList<nBits>& faceValues,
                const CombineOp& cop
            );

            template <unsigned nBits>
            static void swapFaceList
            (
                const polyMesh& mesh,
                PackedList<nBits>& faceValues
            );

            template <unsigned nBits, class CombineOp>
            static void syncPointList
            (
                const polyMesh& mesh,
                PackedList<nBits>& pointValues,
                const CombineOp& cop,
                const unsigned int nullValue
            );

            template <unsigned nBits, class CombineOp>
            static void syncEdgeList
            (
                const polyMesh& mesh,
                PackedList<nBits>& edgeValues,
                const CombineOp& cop,
                const unsigned int nullValue
            );

        // Other

            //- Get per point whether is it master (of a coupled set of points)
            static PackedBoolList getMasterPoints(const polyMesh&);

            //- Get per edge whether is it master (of a coupled set of edges)
            static PackedBoolList getMasterEdges(const polyMesh&);

            //- Get per face whether is it master (of a coupled set of faces)
            static PackedBoolList getMasterFaces(const polyMesh&);

};


template<>
void syncTools::transform::operator()
(
    const coupledPolyPatch&,
    Field<label>&
) const;
template<>
void syncTools::transform::operator()
(
    const coupledPolyPatch&,
    Map<label>&
) const;
template<>
void syncTools::transform::operator()
(
    const coupledPolyPatch&,
    EdgeMap<label>&
) const;

template<>
void syncTools::transform::operator()
(
    const coupledPolyPatch& cpp,
    Field<bool>& fld
) const;
template<>
void syncTools::transform::operator()
(
    const coupledPolyPatch&,
    Map<bool>&
) const;
template<>
void syncTools::transform::operator()
(
    const coupledPolyPatch&,
    EdgeMap<bool>&
) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "syncToolsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
