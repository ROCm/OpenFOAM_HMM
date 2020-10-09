/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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
    Foam::treeBoundBox

Description
    Standard boundBox with extra functionality for use in octree.

    Numbering of corner points is according to octant numbering.

    On the back plane (z=0):

    \verbatim
        Y
        ^
        |
        +--------+
        |2      3|
        |        |
        |        |
        |        |
        |0      1|
        +--------+->X
    \endverbatim

    For the front plane add 4 to the point labels.

Note
    When a bounding box is created without any points, it creates an inverted
    bounding box. Points can be added later and the bounding box will grow to
    include them.

SourceFiles
    treeBoundBoxI.H
    treeBoundBox.C
    treeBoundBoxTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef treeBoundBox_H
#define treeBoundBox_H

#include "boundBox.H"
#include "direction.H"
#include "pointField.H"
#include "faceList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Random;
class treeBoundBox;
Istream& operator>>(Istream& is, treeBoundBox& bb);
Ostream& operator<<(Ostream& os, const treeBoundBox& bb);

/*---------------------------------------------------------------------------*\
                        Class treeBoundBox Declaration
\*---------------------------------------------------------------------------*/

class treeBoundBox
:
    public boundBox
{
public:

    // Enumerations

        //- Bits used for octant/point encoding.
        //  Every octant/corner point is the combination of three directions.
        enum octantBit : direction
        {
            RIGHTHALF = 0x1,  //!< 1: positive x-direction
            TOPHALF   = 0x2,  //!< 2: positive y-direction
            FRONTHALF = 0x4   //!< 4: positive z-direction
        };

        //- Face codes. Identical order and meaning as per hex cellmodel
        //- and boundBox, but have a different point ordering.
        enum faceId
        {
            LEFT   = 0,  //!< 0: x-min, left
            RIGHT  = 1,  //!< 1: x-max, right
            BOTTOM = 2,  //!< 2: y-min, bottom
            TOP    = 3,  //!< 3: y-max, top
            BACK   = 4,  //!< 4: z-min, back
            FRONT  = 5   //!< 5: z-max, front
        };

        //- Bits used for face encoding
        enum faceBit : direction
        {
            NOFACE    = 0,
            LEFTBIT   = 0x1 << LEFT,    //!<  1: x-min, left
            RIGHTBIT  = 0x1 << RIGHT,   //!<  2: x-max, right
            BOTTOMBIT = 0x1 << BOTTOM,  //!<  4: y-min, bottom
            TOPBIT    = 0x1 << TOP,     //!<  8: y-max, top
            BACKBIT   = 0x1 << BACK,    //!< 16: z-min, back
            FRONTBIT  = 0x1 << FRONT    //!< 32: z-max, front
        };

        //- Edges codes.
        //  E01 = edge between 0 and 1.
        enum edgeId
        {
            E01 = 0,
            E13 = 1,
            E23 = 2,
            E02 = 3,

            E45 = 4,
            E57 = 5,
            E67 = 6,
            E46 = 7,

            E04 = 8,
            E15 = 9,
            E37 = 10,
            E26 = 11
        };


    // Static Data Members

        //- Face to point addressing
        static const faceList faces;

        //- Edge to point addressing
        static const edgeList edges;


    // Constructors

        //- Construct without any points - an inverted bounding box
        inline treeBoundBox();

        //- Construct from a boundBox
        inline explicit treeBoundBox(const boundBox& bb);

        //- Construct a bounding box containing a single initial point
        inline explicit treeBoundBox(const point& pt);

        //- Construct from components
        inline treeBoundBox(const point& min, const point& max);

        //- Construct as the bounding box of the given pointField.
        //  Local processor domain only (no reduce as in boundBox)
        explicit treeBoundBox(const UList<point>& points);

        //- Construct as subset of points
        //  Local processor domain only (no reduce as in boundBox)
        treeBoundBox(const UList<point>& points, const labelUList& indices);

        //- Construct as subset of points
        //  The indices could be from edge/triFace etc.
        //  Local processor domain only (no reduce as in boundBox)
        template<unsigned N>
        treeBoundBox
        (
            const UList<point>& points,
            const FixedList<label, N>& indices
        );


        //- Construct from Istream
        inline treeBoundBox(Istream& is);


    // Member functions

        // Access

            //- Typical dimension length,height,width
            inline scalar typDim() const;

            //- Vertex coordinates. In octant coding.
            tmp<pointField> points() const;


        // Check

            //- Corner point of given octant
            inline point corner(const direction octant) const;

            //- Sub-box of given octant. Midpoint calculated.
            treeBoundBox subBbox(const direction octant) const;

            //- Sub-box given by octant number. Midpoint provided.
            treeBoundBox subBbox(const point& mid, const direction) const;

            //- Returns octant number given point and the calculated midpoint.
            inline direction subOctant
            (
                const point& pt
            ) const;

            //- Returns octant number given point and midpoint.
            static inline direction subOctant
            (
                const point& mid,
                const point& pt
            );

            //- Returns octant number given point and the calculated midpoint.
            //  onEdge set if the point is on edge of subOctant
            inline direction subOctant
            (
                const point& pt,
                bool& onEdge
            ) const;

            //- Returns octant number given point and midpoint.
            //  onEdge set if the point is on edge of subOctant
            static inline direction subOctant
            (
                const point& mid,
                const point& pt,
                bool& onEdge
            );

            //- Returns octant number given intersection and midpoint.
            //  onEdge set if the point is on edge of subOctant
            //  If onEdge, the direction vector determines which octant to use
            //  (acc. to which octant the point would be if it were moved
            //  along dir)
            static inline direction subOctant
            (
                const point& mid,
                const vector& dir,
                const point& pt,
                bool& onEdge
            );

            //- Calculates optimal order to look for nearest to point.
            //  First will be the octant containing the point,
            //  second the octant with boundary nearest to the point etc.
            inline void searchOrder
            (
                 const point& pt,
                 FixedList<direction, 8>& octantOrder
            ) const;

            //- Overlaps other bounding box?
            using boundBox::overlaps;

            //- Intersects segment; set point to intersection position and face,
            //  return true if intersection found.
            //  (pt argument used during calculation even if not intersecting).
            //  Calculates intersections from outside supplied vector
            //  (overallStart, overallVec). This is so when
            //  e.g. tracking through lots of consecutive boxes
            // (typical octree) we're not accumulating truncation errors. Set
            // to start, (end-start) if not used.
            bool intersects
            (
                const point& overallStart,
                const vector& overallVec,
                const point& start,
                const point& end,
                point& pt,
                direction& ptBits
            ) const;

            //- Like above but does not return faces point is on
            bool intersects
            (
                const point& start,
                const point& end,
                point& pt
            ) const;

            //- Contains point or other bounding box?
            using boundBox::contains;

            //- Contains point (inside or on edge) and moving in direction
            //  dir would cause it to go inside.
            bool contains(const vector& dir, const point&) const;

            //- Code position of point on bounding box faces
            direction faceBits(const point& pt) const;

            //- Position of point relative to bounding box
            direction posBits(const point& pt) const;

            //- Calculate nearest and furthest (to point) vertex coords of
            //  bounding box
            void calcExtremities
            (
                const point& pt,
                point& nearest,
                point& furthest
            ) const;

            //- Returns distance point to furthest away corner.
            scalar maxDist(const point& pt) const;

            //- Compare distance to point with other bounding box
            //  return:
            //  -1 : all vertices of my bounding box are nearer than any of
            //       other
            //  +1 : all vertices of my bounding box are further away than
            //       any of other
            //   0 : none of the above.
            label distanceCmp(const point& pt, const treeBoundBox& other) const;

            //- Return slightly wider bounding box
            //  Extends all dimensions with s*span*Random::sample01<scalar>()
            //  and guarantees in any direction s*mag(span) minimum width
            inline treeBoundBox extend(Random& rndGen, const scalar s) const;


    // IOstream operator

        friend Istream& operator>>(Istream& is, treeBoundBox& bb);
        friend Ostream& operator<<(Ostream& os, const treeBoundBox& bb);
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Contiguous data for treeBoundBox
template<> struct is_contiguous<treeBoundBox> : is_contiguous<boundBox> {};

//- Contiguous scalar data for treeBoundBox
template<> struct is_contiguous_scalar<treeBoundBox>
:
    is_contiguous_scalar<boundBox>
{};


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

inline bool operator==(const treeBoundBox& a, const treeBoundBox& b);
inline bool operator!=(const treeBoundBox& a, const treeBoundBox& b);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "treeBoundBoxI.H"

#ifdef NoRepository
    #include "treeBoundBoxTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
