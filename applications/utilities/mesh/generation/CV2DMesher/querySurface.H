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
    querySurface

Description
    Searchable triSurface using an octree to speed-up queries.

SourceFiles
    querySurface.C

\*---------------------------------------------------------------------------*/

#ifndef querySurface_H
#define querySurface_H

#include "triSurface.H"
#include "treeDataTriSurface.H"
#include "indexedOctree.H"
#include "Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class querySurface Declaration
\*---------------------------------------------------------------------------*/

class querySurface
:
    public triSurface
{
    // Private data

        Random rndGen_;

        // Bounding box of surface. Used for relative tolerances.
        treeBoundBox bb_;

        // Search engine on surface
        indexedOctree<treeDataTriSurface> tree_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        querySurface(const querySurface&);

        //- Disallow default bitwise assignment
        void operator=(const querySurface&);


public:

    // Constructors

        //- Construct given file name of the surface
        querySurface(const fileName& surfaceFileName);


    // Destructor

        ~querySurface();


    // Member Functions

        // Access

            const treeBoundBox& bb() const
            {
                return bb_;
            }

            const indexedOctree<treeDataTriSurface>& tree() const
            {
                return tree_;
            }


        // Query

            //- Extract feature edges/points(2D)
            //  using the given feature angle in deg
            labelList extractFeatures2D(const scalar featAngle) const;

            //- Returns inside, outside or mixed (= close to surface)
            indexedOctree<Foam::treeDataTriSurface>::volumeType insideOutside
            (
                const scalar searchSpan2,
                const point& pt
            ) const;

            //- Check if point is inside surface
            bool inside(const point& pt) const;

            //- Check if point is outside surface
            bool outside(const point& pt) const;

            //- Check if point is inside surface by at least dist2
            bool wellInside(const point& pt, const scalar dist2) const;

            //- Check if point is outside surface by at least dist2
            bool wellOutside(const point& pt, const scalar dist2) const;


        // Write

            void writeTreeOBJ() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "querySurfaceI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
