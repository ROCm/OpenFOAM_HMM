/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
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
    A bounding box defined in terms of the points at its extremities.

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

// Forward declaration of friend functions and operators

class boundBox;
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

    // Static data members

        //- A large boundBox: min/max == -/+ ROOTVGREAT
        static const boundBox greatBox;

        //- A large inverted boundBox: min/max == +/- ROOTVGREAT
        static const boundBox invertedBox;

        //- Faces to point addressing, as per a 'hex' cell
        static const faceList faces;


    // Constructors

        //- Construct without any points - an inverted bounding box
        inline boundBox();

        //- Construct a bounding box containing a single initial point
        explicit inline boundBox(const point& pt);

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
        template<unsigned Size>
        boundBox
        (
            const UList<point>& points,
            const FixedList<label, Size>& indices,
            bool doReduce = true
        );

        //- Construct from Istream
        inline boundBox(Istream& is);


    // Member functions

        // Access

            //- Bounding box is inverted, contains no points.
            inline bool empty() const;

            //- Clear bounding box of all points - make it an inverted box
            inline void clear();

            //- Minimum describing the bounding box
            inline const point& min() const;

            //- Maximum describing the bounding box
            inline const point& max() const;

            //- Minimum describing the bounding box, non-const access
            inline point& min();

            //- Maximum describing the bounding box, non-const access
            inline point& max();

            //- The midpoint of the bounding box
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

            //- Corner points in an order corresponding to a 'hex' cell
            tmp<pointField> points() const;


        // Manipulate

            //- Extend to include the second box.
            inline void add(const boundBox& bb);

            //- Extend to include the point.
            inline void add(const point& pt);

            //- Extend to include the points.
            inline void add(const UList<point>& points);

            //- Extend to include the points from the temporary point field.
            inline void add(const tmp<pointField>& tpoints);

            //- Extend to include the subset of the point field.
            //  The indices could be from cell/face etc.
            inline void add
            (
                const UList<point>& points,
                const labelUList& indices
            );

            //- Extend to include the points.
            template<unsigned Size>
            void add(const FixedList<point, Size>& points);

            //- Extend to include the subset of the point field.
            //  The indices could be from edge/triFace etc.
            template<unsigned Size>
            void add
            (
                const UList<point>& points,
                const FixedList<label, Size>& indices
            );

            //- Inflate box by factor*mag(span) in all dimensions
            void inflate(const scalar s);

            //- Parallel reduction of min/max values
            void reduce();

            //- Intersection (union) with the second box.
            //  The return value is true if the intersection is non-empty.
            bool intersect(const boundBox& bb);


        // Query

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

            //- Contains all of the points? (inside or on edge)
            bool contains(const UList<point>& points) const;

            //- Contains all of the subset of points? (inside or on edge)
            bool contains
            (
                const UList<point>& points,
                const labelUList& indices
            ) const;

            //- Contains all of the subset of points? (inside or on edge)
            template<unsigned Size>
            bool contains
            (
                const UList<point>& points,
                const FixedList<label, Size>& indices
            ) const;


            //- Contains any of the points? (inside or on edge)
            bool containsAny(const UList<point>& points) const;

            //- Contains any of the subset of points? (inside or on edge)
            bool containsAny
            (
                const UList<point>& points,
                const labelUList& indices
            ) const;

            //- Contains any of the subset of points? (inside or on edge)
            template<unsigned Size>
            bool containsAny
            (
                const UList<point>& points,
                const FixedList<label, Size>& indices
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


//- Data associated with boundBox type are contiguous
template<>
inline bool contiguous<boundBox>() {return contiguous<point>();}

// Global Operators

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
