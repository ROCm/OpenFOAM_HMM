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
    Foam::triSurfaceMesh

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
    public objectRegistry,      // so we can store fields
    public triSurface
{
private:

    // Private member data

        //- Search tree (triangles)
        mutable autoPtr<indexedOctree<treeDataTriSurface> > tree_;

        //- Search tree for boundary edges.
        mutable autoPtr<indexedOctree<treeDataEdge> > edgeTree_;

        //- Names of regions
        mutable wordList regions_;

        //- Is surface closed
        mutable label surfaceClosed_;

    // Private Member Functions

        //- Check file existence
        static const fileName& checkFile
        (
            const fileName& fName,
            const fileName& objectName
        );

        //- Check whether surface is closed without calculating any permanent
        //  addressing.
        bool isSurfaceClosed() const;

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

        //- Construct from dictionary (used by searchableSurface)
        triSurfaceMesh
        (
            const IOobject& io,
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

            virtual const wordList& regions() const;

            //- Whether supports volume type below. I.e. whether is closed.
            virtual bool hasVolumeType() const;


        // Multiple point queries.

            virtual void findNearest
            (
                const pointField& sample,
                const scalarField& nearestDistSqr,
                List<pointIndexHit>&
            ) const;

            virtual void findLine
            (
                const pointField& start,
                const pointField& end,
                List<pointIndexHit>&
            ) const;

            virtual void findLineAny
            (
                const pointField& start,
                const pointField& end,
                List<pointIndexHit>&
            ) const;

            //- Get all intersections in order from start to end.
            virtual void findLineAll
            (
                const pointField& start,
                const pointField& end,
                List<List<pointIndexHit> >&
            ) const;

            //- From a set of points and indices get the region
            virtual void getRegion
            (
                const List<pointIndexHit>&,
                labelList& region
            ) const;

            //- From a set of points and indices get the normal
            virtual void getNormal
            (
                const List<pointIndexHit>&,
                vectorField& normal
            ) const;

            //- Determine type (inside/outside/mixed) for point. unknown if
            //  cannot be determined (e.g. non-manifold surface)
            virtual void getVolumeType
            (
                const pointField&,
                List<volumeType>&
            ) const;


        // regIOobject implementation

            bool writeData(Ostream&) const
            {
                notImplemented("triSurfaceMesh::writeData(Ostream&) const");
                return false;
            }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
