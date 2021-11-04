/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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
    Foam::line

Description
    A line primitive.

SourceFiles
    lineI.H

\*---------------------------------------------------------------------------*/

#ifndef line_H
#define line_H

#include "vector.H"
#include "PointHit.H"
#include "point2D.H"
#include "FixedList.H"
#include "UList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations

template<class Point, class PointRef> class line;

template<class Point, class PointRef>
inline Istream& operator>>(Istream& is, line<Point, PointRef>& l);

template<class Point, class PointRef>
inline Ostream& operator<<(Ostream& os, const line<Point, PointRef>& l);


/*---------------------------------------------------------------------------*\
                           Class line Declaration
\*---------------------------------------------------------------------------*/

template<class Point, class PointRef>
class line
{
    // Private Data

        //- First point
        PointRef a_;

        //- Second point
        PointRef b_;


public:

    // Constructors

        //- Construct from two points
        inline line(const Point& from, const Point& to);

        //- Construct from two points in the list of points
        //  The indices could be from edge etc.
        inline line
        (
            const UList<Point>& points,
            const FixedList<label, 2>& indices
        );

        //- Construct from Istream
        inline explicit line(Istream& is);


    // Member Functions

    // Access

        //- Return first point
        inline PointRef first() const noexcept;

        //- Return second (last) point
        inline PointRef second() const noexcept;

        //- Return last (second) point
        inline PointRef last() const noexcept;

        //- Return first point
        inline PointRef start() const noexcept;

        //- Return second (last) point
        inline PointRef end() const noexcept;


    // Properties

        //- Return centre (centroid)
        inline Point centre() const;

        //- Return scalar magnitude
        inline scalar mag() const;

        //- Return start-to-end vector
        inline Point vec() const;

        //- Return the unit vector (start-to-end)
        inline Point unitVec() const;

        //- Return nearest distance to line from a given point
        //  If the nearest point is on the line, return a hit
        PointHit<Point> nearestDist(const Point& p) const;

        //- Return nearest distance from line to line. Returns distance
        //  and sets both points (one on *this, one on the provided
        //  linePointRef.
        scalar nearestDist
        (
            const line<Point, const Point&>& edge,
            Point& thisPoint,
            Point& edgePoint
        ) const;


    // IOstream operators

        friend Istream& operator>> <Point, PointRef>
        (
            Istream& is,
            line& l
        );

        friend Ostream& operator<< <Point, PointRef>
        (
            Ostream& os,
            const line& l
        );
};


//- 2D specialisation
template<>
scalar line<point2D, const point2D&>::nearestDist
(
    const line<point2D, const point2D&>& edge,
    point2D& thisPoint,
    point2D& edgePoint
) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "lineI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
