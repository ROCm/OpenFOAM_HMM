/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2019 OpenCFD Ltd.
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
    Foam::boundBox

Description
    A bounding box defined in terms of min/max extrema points.

Note
    When a bounding box is created without any points, it creates an inverted
    bounding box. Points can be added later and the bounding box will grow to
    include them.

\*---------------------------------------------------------------------------*/

#ifndef boundBox_H
#define boundBox_H

#include "pointField.H"
#include "faceList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class boundBox;
class plane;
template<class T> class tmp;

Istream& operator>>(Istream& is, boundBox& bb);
Ostream& operator<<(Ostream& os, const boundBox& bb);


/*---------------------------------------------------------------------------*\
                          Class boundBox Declaration
\*---------------------------------------------------------------------------*/

class boundBox
{
    // Private data

        //- Minimum and maximum points describing the bounding box
        point min_, max_;

public:

    // Static Data Members

        //- Bits used for (x/y/z) direction encoding.
        enum directionBit : direction
        {
            XDIR = 0x1,  //!< 1: x-direction (vector component 0)
            YDIR = 0x2,  //!< 2: y-direction (vector component 1)
            ZDIR = 0x4   //!< 4: z-direction (vector component 2)
        };

        //- A large boundBox: min/max == -/+ ROOTVGREAT
        static const boundBox greatBox;

        //- A large inverted boundBox: min/max == +/- ROOTVGREAT
        static const boundBox invertedBox;

        //- Faces to point addressing, as per a 'hex' cell
        static const faceList faces;

        //- The unit normal per face
        static const FixedList<vector, 6> faceNormals;


    // Constructors

        //- Construct without any points - an inverted bounding box
        inline boundBox();

        //- Construct a bounding box containing a single initial point
        inline explicit boundBox(const point& pt);

        //- Construct from components
        inline boundBox(const point& min, const point& max);

        //- Construct as the bounding box of the given points
        //  Does parallel communication (doReduce = true)
        explicit boundBox(const UList<point>& points, bool doReduce = true);

        //- Construct as the bounding box of the given temporary pointField.
        //  Does parallel communication (doReduce = true)
        explicit boundBox(const tmp<pointField>& tpoints, bool doReduce = true);

        //- Construct bounding box as an indirect subset of the points.
        //  The indices could be from cell/face etc.
        //  Does parallel communication (doReduce = true)
        boundBox
        (
            const UList<point>& points,
            const labelUList& indices,
            bool doReduce = true
        );

        //- Construct bounding box as an indirect subset of the points.
        //  The indices could be from edge/triFace etc.
        //  Does parallel communication (doReduce = true)
        template<unsigned N>
        boundBox
        (
            const UList<point>& points,
            const FixedList<label, N>& indices,
            bool doReduce = true
        );

        //- Construct from Istream
        inline boundBox(Istream& is);


    // Member Functions

        // Access

            //- Bounding box is inverted, contains no points.
            inline bool empty() const;

            //- Bounding box is non-inverted.
            inline bool valid() const;

            //- Minimum describing the bounding box
            inline const point& min() const;

            //- Maximum describing the bounding box
            inline const point& max() const;

            //- Minimum describing the bounding box, non-const access
            inline point& min();

            //- Maximum describing the bounding box, non-const access
            inline point& max();

            //- The centre (midpoint) of the bounding box
            inline point centre() const;

            //- The midpoint (centre) of the bounding box. Identical to centre()
            inline point midpoint() const;

            //- The bounding box span (from minimum to maximum)
            inline vector span() const;

            //- The magnitude of the bounding box span
            inline scalar mag() const;

            //- The volume of the bound box
            inline scalar volume() const;

            //- Smallest length/height/width dimension
            inline scalar minDim() const;

            //- Largest length/height/width dimension
            inline scalar maxDim() const;

            //- Average length/height/width dimension
            inline scalar avgDim() const;

            //- Count the number of positive, non-zero dimensions.
            //  \return -1 if any dimensions are negative,
            //  0 = 0D (point),
            //  1 = 1D (line aligned with an axis),
            //  2 = 2D (plane aligned with an axis),
            //  3 = 3D (box)
            inline label nDim() const;

            //- Corner points in an order corresponding to a 'hex' cell
            tmp<pointField> points() const;

            //- Face midpoints
            tmp<pointField> faceCentres() const;

            //- Face centre of given face index
            point faceCentre(const direction facei) const;


        // Manipulate

            //- Clear bounding box and make it an inverted box
            inline void clear();

            //- Extend to include the second box.
            inline void add(const boundBox& bb);

            //- Extend to include the point.
            inline void add(const point& pt);

            //- Extend to include the points.
            inline void add(const UList<point>& points);

            //- Extend to include the points from the temporary point field.
            inline void add(const tmp<pointField>& tpoints);

            //- Extend to include the points.
            template<unsigned N>
            void add(const FixedList<point, N>& points);

            //- Extend to include a (subsetted) point field.
            //  The indices could be from edge/triFace etc.
            template<unsigned N>
            void add
            (
                const UList<point>& points,
                const FixedList<label, N>& indices
            );

            //- Extend to include a (subsetted) point field.
            //
            //  \tparam IntContainer  A container with an iterator that
            //      dereferences to an label
            template<class IntContainer>
            void add
            (
                const UList<point>& points,
                const IntContainer& indices
            );

            //- Inflate box by factor*mag(span) in all dimensions
            void inflate(const scalar s);

            //- Parallel reduction of min/max values
            void reduce();


        // Query

            //- Intersection (union) with the second box.
            //  The return value is true if the intersection is non-empty.
            bool intersect(const boundBox& bb);

            //- Does plane intersect this bounding box.
            //  There is an intersection if the plane segments the corner points
            //  \note the results are unreliable when plane coincides almost
            //  exactly with a box face
            bool intersects(const plane& pln) const;

            //- Overlaps/touches boundingBox?
            inline bool overlaps(const boundBox& bb) const;

            //- Overlaps boundingSphere (centre + sqr(radius))?
            inline bool overlaps
            (
                const point& centre,
                const scalar radiusSqr
            ) const;

            //- Contains point? (inside or on edge)
            inline bool contains(const point& pt) const;

            //- Fully contains other boundingBox?
            inline bool contains(const boundBox& bb) const;

            //- Contains point? (inside only)
            inline bool containsInside(const point& pt) const;

            //- Contains all points? (inside or on edge)
            bool contains(const UList<point>& points) const;

            //- Contains all of the (subsetted) points? (inside or on edge)
            template<unsigned N>
            bool contains
            (
                const UList<point>& points,
                const FixedList<label, N>& indices
            ) const;


            //- Contains all of the (subsetted) points? (inside or on edge)
            //
            //  \tparam IntContainer  A container with an iterator that
            //      dereferences to an label
            template<class IntContainer>
            bool contains
            (
                const UList<point>& points,
                const IntContainer& indices
            ) const;


            //- Contains any of the points? (inside or on edge)
            bool containsAny(const UList<point>& points) const;

            //- Contains any of the (subsetted) points? (inside or on edge)
            template<unsigned N>
            bool containsAny
            (
                const UList<point>& points,
                const FixedList<label, N>& indices
            ) const;

            //- Contains any of the (subsetted) points? (inside or on edge)
            //
            //  \tparam IntContainer  A container with an iterator that
            //      dereferences to an label
            template<class IntContainer>
            bool containsAny
            (
                const UList<point>& points,
                const IntContainer& indices
            ) const;


            //- Return the nearest point on the boundBox to the supplied point.
            //  If point is inside the boundBox then the point is returned
            //  unchanged.
            point nearest(const point& pt) const;


    // Member Operators

        //- Extend box to include the second box, as per the add() method.
        inline void operator+=(const boundBox& bb);


    // IOstream operator

        friend Istream& operator>>(Istream& is, boundBox& bb);
        friend Ostream& operator<<(Ostream& os, const boundBox& bb);
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Contiguous data for boundBox
template<> struct is_contiguous<boundBox> : is_contiguous<point> {};

//- Contiguous scalar data for boundBox
template<> struct is_contiguous_scalar<boundBox>
:
    is_contiguous_scalar<point>
{};


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

inline bool operator==(const boundBox& a, const boundBox& b);
inline bool operator!=(const boundBox& a, const boundBox& b);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "boundBoxI.H"

#ifdef NoRepository
    #include "boundBoxTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
