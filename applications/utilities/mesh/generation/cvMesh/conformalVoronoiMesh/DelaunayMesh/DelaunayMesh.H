/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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
    Foam::DelaunayMesh

Description
    The vertex and cell classes must have an index defined

SourceFiles
    DelaunayMeshI.H
    DelaunayMesh.C
    DelaunayMeshIO.C

\*---------------------------------------------------------------------------*/

#ifndef DelaunayMesh_H
#define DelaunayMesh_H

#include "Pair.H"
#include "HashSet.H"
#include "FixedList.H"
#include "boundBox.H"
#include "indexedVertex.H"
#include "CGALTriangulation3Ddefs.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                        Class DelaunayMesh Declaration
\*---------------------------------------------------------------------------*/

template<class Triangulation>
class DelaunayMesh
:
    public Triangulation
{
public:

        typedef typename Triangulation::Cell_handle Cell_handle;
        typedef typename Triangulation::Vertex_handle Vertex_handle;
        typedef typename Triangulation::Point Point;
        typedef typename Triangulation::Facet Facet;

        typedef typename Triangulation::Finite_vertices_iterator
            Finite_vertices_iterator;
        typedef typename Triangulation::Finite_cells_iterator
            Finite_cells_iterator;
        typedef typename Triangulation::Finite_facets_iterator
            Finite_facets_iterator;

        typedef HashSet
        <
           Pair<label>,
           FixedList<label, 2>::Hash<>
        > labelPairHashSet;


private:

    // Private data

        //- Keep track of the number of vertices that have been added.
        //  This allows a unique index to be assigned to each vertex.
        mutable label vertexCount_;

        //- Keep track of the number of cells that have been added.
        //  This allows a unique index to be assigned to each cell.
        mutable label cellCount_;

        //- Spatial sort traits to use with a pair of point pointers and an int.
        //  Taken from a post on the CGAL lists: 2010-01/msg00004.html by
        //  Sebastien Loriot (Geometry Factory).
        struct Traits_for_spatial_sort
        :
            public Triangulation::Geom_traits
        {
            typedef typename Triangulation::Geom_traits Gt;

            typedef std::pair<const typename Triangulation::Point*, int>
                Point_3;

            struct Less_x_3
            {
                bool operator()(const Point_3& p, const Point_3& q) const;
            };

            struct Less_y_3
            {
                bool operator()(const Point_3& p, const Point_3& q) const;
            };

            struct Less_z_3
            {
                bool operator()(const Point_3& p, const Point_3& q) const;
            };

            Less_x_3 less_x_3_object() const;
            Less_y_3 less_y_3_object() const;
            Less_z_3 less_z_3_object() const;
        };


    // Private Member Functions

        void sortFaces
        (
            faceList& faces,
            labelList& owner,
            labelList& neighbour
        ) const;

        void addPatches
        (
            const label nInternalFaces,
            faceList& faces,
            labelList& owner,
            PtrList<dictionary>& patchDicts,
            const List<DynamicList<face> >& patchFaces,
            const List<DynamicList<label> >& patchOwners
        ) const;

        //- Disallow default bitwise copy construct
        DelaunayMesh(const DelaunayMesh<Triangulation>&);

        //- Disallow default bitwise assignment
        void operator=(const DelaunayMesh<Triangulation>&);


public:

    // Constructors

        //- Construct from components
        DelaunayMesh();


    //- Destructor
    ~DelaunayMesh();


    // Member Functions

        inline label getNewVertexIndex() const;

        inline label getNewCellIndex() const;

        inline label cellCount() const;

        inline void resetCellCount();

        inline label vertexCount() const;

        inline void resetVertexCount();


        //- Remove the entire triangulation
        void reset();

        void insertPoints(const List<Vb>& vertices);

        //- Function inserting points into a triangulation and setting the
        //  index and type data of the point in the correct order. This is
        //  faster than inserting points individually.
        //
        //  Adapted from a post on the CGAL lists: 2010-01/msg00004.html by
        //  Sebastien Loriot (Geometry Factory).
        template<class PointIterator>
        void rangeInsertWithInfo
        (
            PointIterator begin,
            PointIterator end,
            bool printErrors = true
        );


    // Queries

        void printInfo(Ostream& os) const;

        void printVertexInfo(Ostream& os) const;

        //- Create an fvMesh from the triangulation.
        //  The mesh is not parallel consistent - only used for viewing
        autoPtr<fvMesh> createMesh
        (
            const fileName& name,
            const Time& runTime,
            labelList& vertexMap,
            labelList& cellMap
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DelaunayMeshI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "DelaunayMesh.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
