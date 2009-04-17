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
    Foam::featureEdgeMesh

Description

    Feature points are a sorted subset at the start of the overall points list:
        0 .. concaveStart_-1                : convex points (w.r.t normals)
        concaveStart_-1 .. mixedStart_-1    : concave points
        mixedStart_ .. nonFeatureStart_     : mixed internal/external points
        nonFeatureStart_ .. size-1          : non-feature points

    Feature edges are the edgeList of the edgeMesh and are sorted:
        0 .. internalStart_-1           : external edges (convex)
        internalStart_ .. flatStart_-1  : internal edges (concave)
        flatStart_ .. openStart_-1      : flat edges (neither concave or convex)
                                          can arise from region interfaces on
                                          flat surfaces
        openStart_ .. multipleStart_-1  : open edges (e.g. from baffle surfaces)
        multipleStart_ .. size-1        : multiply connected edges


SourceFiles
    featureEdgeMeshI.H
    featureEdgeMesh.C

\*---------------------------------------------------------------------------*/

#ifndef featureEdgeMesh_H
#define featureEdgeMesh_H

#include "edgeMesh.H"
#include "indexedOctree.H"
#include "treeDataEdge.H"
#include "pointIndexHit.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class featureEdgeMesh Declaration
\*---------------------------------------------------------------------------*/

class featureEdgeMesh
:
    public edgeMesh
{

public:

        enum pointStatus
        {
            CONVEX,     // Fully convex point (w.r.t normals)
            CONCAVE,    // Fully concave point
            MIXED,      // A point surrounded by both convex and concave edges
            NONFEATURE  // Not a feature point
        };

        enum edgeStatus
        {
            EXTERNAL,   // "Convex" edge
            INTERNAL,   // "Concave" edge
            FLAT,       // Neither concave or convex, on a flat surface
            OPEN,       // i.e. only connected to one face
            MULTIPLE,   // Multiply connected (connected to more than two faces)
            NONE        // Not a feature edge (consistency with surfaceFeatures)
        };

private:

    // Private data

        //- Index of the start of the concave feature points
        label concaveStart_;

        //- Index of the start of the mixed type feature points
        label mixedStart_;

        //- Index of the start of the non-feature points
        label nonFeatureStart_;

        //- Index of the start of the internal feature edges
        label internalStart_;

        //- Index of the start of the flat feature edges
        label flatStart_;

        //- Index of the start of the open feature edges
        label openStart_;

        //- Index of the start of the multiply-connected feature edges
        label multipleStart_;

        //- Normals of the features, to be referred to by index by both feature
        //  points and edges, unsorted
        vectorField normals_;

        //- Indices of the normals that are adjacent to the feature edges
        labelListList edgeNormals_;

        //- Indices of the normals that are adjacent to the feature points
        labelListList featurePointNormals_;

        //- Indices of all edges, pre-calculated and stored
        labelField allEdges_;

        //- Search tree for edges.
        mutable autoPtr<indexedOctree<treeDataEdge> > edgeTree_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        featureEdgeMesh(const featureEdgeMesh&);

        //- Disallow default bitwise assignment
        void operator=(const featureEdgeMesh&);


public:

    // Constructors

        //- Construct (read) given an IOobject
        featureEdgeMesh(const IOobject&);

        //- Construct from edgeMesh data and
        featureEdgeMesh
        (
            const IOobject& io,
            const pointField& pts,
            const edgeList& eds,
            const vectorField& normals
        );


    //- Destructor
    ~featureEdgeMesh();


    // Member Functions

        // Find

            //- Find nearest surface edge for each sample point.
            void nearestFeatureEdge
            (
                const pointField& samples,
                scalarField searchDistSqr,
                labelList& edgeLabel,
                pointField& edgePoint,
                labelListList& adjacentNormals
            ) const;

        // Access

            //- Return the normals of the surfaces adjacent to the feature edges
            //  and points
            inline const vectorField& normals() const;

            //- Demand driven contruction of octree for boundary edges
            const indexedOctree<treeDataEdge>& edgeTree() const;

        // Check

        // Edit

        // Write
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "featureEdgeMeshI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
