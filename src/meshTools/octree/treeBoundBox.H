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
    Foam::treeBoundBox

Description
    Standard boundBox with extra functionality for use in octree.

SourceFiles
    treeBoundBoxI.H
    treeBoundBox.C

\*---------------------------------------------------------------------------*/

#ifndef treeBoundBox_H
#define treeBoundBox_H

#include "boundBox.H"
#include "direction.H"
#include "pointField.H"
#include "edgeList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Random;

/*---------------------------------------------------------------------------*\
                           Class treeBoundBox Declaration
\*---------------------------------------------------------------------------*/

class treeBoundBox
:
    public boundBox
{

public:

    // Static data members

        static treeBoundBox greatBox;

        // Bits used for octant coding
        static const direction RIGHTHALF   = 0x1 << 0;
        static const direction TOPHALF     = 0x1 << 1;
        static const direction FRONTHALF   = 0x1 << 2;

        // Bits used for face coding
        static const direction NOFACE     = 0;
        static const direction LEFTBIT    = 0x1 << 0;
        static const direction RIGHTBIT   = 0x1 << 1;
        static const direction BELOWBIT   = 0x1 << 2;
        static const direction ABOVEBIT   = 0x1 << 3;
        static const direction BEHINDBIT  = 0x1 << 4;
        static const direction INFRONTBIT = 0x1 << 5;

        //- Face on which neighbour is
        static direction neighbourFaceBits(const label&);


    // Constructors

        //- Construct null setting points to zero
        inline treeBoundBox();

        //- Construct from components
        inline treeBoundBox(const point& min, const point& max);

        //- Construct from components
        inline treeBoundBox(const boundBox& bb);

        //- Construct as the bounding box of the given pointField. Local
        //  processor domain only (no reduce as in boundBox)
        treeBoundBox(const pointField& points);

        //- Construct from Istream
        treeBoundBox(Istream&);


    // Member functions

        // Access

            //- Smallest of length,height,width
            scalar minDim() const;

            //- Largest of length,height,width
            scalar maxDim() const;

            //- Average of length,height,width
            scalar avgDim() const;

            //- Typical dimension length,height,width
            scalar typDim() const;

            //- vertex coordinates
            pointField points() const;

            //- edges expressed in terms of points()
            edgeList edges() const;


        // Check

            //- Calculates midpoint
            point mid() const;

            //- Sub box given by octant number. Midpoint calculated.
            treeBoundBox subBbox(const direction) const;

            //- Sub box given by octant number. Midpoint provided.
            treeBoundBox subBbox(const point& mid, const direction) const;

            //- Returns octant number given point. Midpoint calculated.
            inline direction subOctant
            (
                const point& sample
            ) const;

            //- Returns octant number given point. Midpoint provided.
            static inline direction subOctant
            (
                const point& mid,
                const point& sample
            );

            //- Returns octant number given point. Midpoint calculated.
            //  onEdge set if sample on edge of subOctant
            inline direction subOctant
            (
                const point& mid,
                bool& onEdge
            ) const;

            //- Returns octant number given point. Midpoint provided.
            //  onEdge set if sample on edge of subOctant
            static inline direction subOctant
            (
                const point& mid,
                const point& sample,
                bool& onEdge
            );

            //- Returns octant number given intersection. Midpoint provided.
            //  onEdge set if sample on edge of subOctant. If onEdge
            // the direction vector determines which octant to use
            // (acc. to which octant the sample would be if it were moved 
            //  along dir)
            static inline direction subOctant
            (
                const point& mid,
                const vector& dir,
                const point& sample,
                bool& onEdge
            );

            //- Intersects other boundingbox?
            inline bool intersects(const treeBoundBox&) const;

            //- Intersects segment; set point to intersection position,
            //  return true if intersection found.
            // (intPt argument used during calculation even if not intersecting)
            bool intersects(const point&, const point&, point& intPt) const;

            //- fully contains bb
            bool contains(const treeBoundBox& bb) const;

            //- Contains point? (inside or on edge)
            inline bool contains(const point&) const;

            //- Contains point? (only inside)
            bool containsNarrow(const point&) const;

            //- Contains point inside or
            //  on edge and moving in direction dir would cause it to go
            //  inside.
            bool contains(const vector& dir, const point&) const;

            //- Position of point relative to bb
            direction posBits(const point&) const;

            //- Calculate nearest and furthest (to sample) vertex coords of
            //  bounding box
            void calcExtremities
            (
                const point& sample,
                point& nearest,
                point& furthest
            ) const;

            //- Returns distance sample to furthest away corner.
            scalar maxDist(const point& sample) const;

            //- Compare distance to point with other bounding box
            //  return:
            //  -1 : all vertices of my bounding box are nearer than any of
            //       other
            //  +1 : all vertices of my bounding box are further away than
            //       any of other
            //   0 : none of the above.
            label distanceCmp(const point&, const treeBoundBox& other) const;

            //- Return slightly wider bounding box
            //  Extends all dimensions with s*span*Random::scalar01()
            inline treeBoundBox extend(Random&, const scalar s) const;

    // Friend Operators

        friend bool operator==(const treeBoundBox&, const treeBoundBox&);
        friend bool operator!=(const treeBoundBox&, const treeBoundBox&);


    // IOstream operator

        friend Istream& operator>>(Istream&, treeBoundBox&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "treeBoundBoxI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
