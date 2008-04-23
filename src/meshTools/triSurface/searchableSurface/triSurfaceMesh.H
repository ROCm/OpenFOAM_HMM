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
    triSurfaceMesh

Description
    IOoject and searching on triSurface

SourceFiles
    triSurfaceMesh.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceMesh_H
#define triSurfaceMesh_H

#include "searchableSurface.H"
#include "objectRegistry.H"
#include "indexedOctree.H"
#include "treeDataTriSurface.H"
#include "treeDataEdge.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class triSurfaceMesh Declaration
\*---------------------------------------------------------------------------*/

class triSurfaceMesh
:
    public searchableSurface,
    public objectRegistry,
    public triSurface
{
private:

    // Private member data

        //- Search tree (triangles)
        mutable autoPtr<indexedOctree<treeDataTriSurface> > tree_;

        //- Search tree for boundary edges.
        mutable autoPtr<indexedOctree<treeDataEdge> > edgeTree_;

    // Private Member Functions

        //- Check file existence
        static const fileName& checkFile
        (
            const fileName& fName,
            const fileName& objectName
        );

        //- Disallow default bitwise copy construct
        triSurfaceMesh(const triSurfaceMesh&);

        //- Disallow default bitwise assignment
        void operator=(const triSurfaceMesh&);


public:

    //- Runtime type information
    TypeName("triSurfaceMesh");


    // Constructors

        //- Construct from triSurface
        triSurfaceMesh(const IOobject&, const triSurface&);

        //- Construct read
        triSurfaceMesh(const IOobject& io);

        //- Construct as searchableSurface
        triSurfaceMesh
        (
            const word& name,
            const objectRegistry& obj,
            const dictionary& dict
        );


    // Destructor

        virtual ~triSurfaceMesh();

        //- Clear storage
        void clearOut();


    // Member Functions

        //- Move points
        virtual void movePoints(const pointField&);

        //- Demand driven contruction of octree
        const indexedOctree<treeDataTriSurface>& tree() const;

        //- Demand driven contruction of octree for boundary edges
        const indexedOctree<treeDataEdge>& edgeTree() const;


        // searchableSurface implementation

            //- Calculate nearest point on surface. Returns
            //  - bool : any point found nearer than nearestDistSqr
            //  - label: relevant index in surface
            //  - point: actual nearest point found
            virtual pointIndexHit findNearest
            (
                const point& sample,
                const scalar nearestDistSqr
            ) const;

            //- Calculate nearest point on edge. Returns
            //  - bool : any point found nearer than nearestDistSqr
            //  - label: relevant index in surface
            //  - point: actual nearest point found
            virtual pointIndexHit findNearestOnEdge
            (
                const point& sample,
                const scalar nearestDistSqr
            ) const;

            //- Find nearest to line. Returns
            //  - bool : any point found?
            //  - label: relevant index in shapes
            //  - point: actual nearest point found
            //  sets:
            //  - tightest  : bounding box
            //  - linePoint : corresponding nearest point on line
            virtual pointIndexHit findNearest
            (
                const linePointRef& ln,
                treeBoundBox& tightest,
                point& linePoint
            ) const;

            //- Find nearest intersection of line between start and end.
            virtual pointIndexHit findLine
            (
                const point& start,
                const point& end
            ) const;

            //- Find any intersection of line between start and end.
            virtual pointIndexHit findLineAny
            (
                const point& start,
                const point& end
            ) const;

            //- Determine type (inside/outside/mixed) for point. unknown if
            //  cannot be determined (e.g. non-manifold surface)
            virtual volumeType getVolumeType(const point&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
