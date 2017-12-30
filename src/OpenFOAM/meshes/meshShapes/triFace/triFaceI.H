/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "IOstreams.H"
#include "face.H"
#include "triPointRef.H"
#include "Swap.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

inline int Foam::triFace::compare(const triFace& a, const triFace& b)
{
    if
    (
        (a[0] == b[0] && a[1] == b[1] && a[2] == b[2])
     || (a[0] == b[1] && a[1] == b[2] && a[2] == b[0])
     || (a[0] == b[2] && a[1] == b[0] && a[2] == b[1])
    )
    {
        // identical
        return 1;
    }
    else if
    (
        (a[0] == b[2] && a[1] == b[1] && a[2] == b[0])
     || (a[0] == b[1] && a[1] == b[0] && a[2] == b[2])
     || (a[0] == b[0] && a[1] == b[2] && a[2] == b[1])
    )
    {
        // same face, but reversed orientation
        return -1;
    }

    return 0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::triFace::triFace()
:
    FixedList<label, 3>(-1)
{}


inline Foam::triFace::triFace
(
    const label a,
    const label b,
    const label c
)
{
    operator[](0) = a;
    operator[](1) = b;
    operator[](2) = c;
}


inline Foam::triFace::triFace(const labelUList& lst)
:
    FixedList<label, 3>(lst)
{}


inline Foam::triFace::triFace(std::initializer_list<label> lst)
:
    FixedList<label, 3>(lst)
{}


inline Foam::triFace::triFace(Istream& is)
:
    FixedList<label, 3>(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::triFace::collapse()
{
    // Cannot resize FixedList, so mark duplicates with '-1'
    // (the lower vertex is retained)
    // catch any '-1' (eg, if called multiple times)

    label n = 3;
    if (operator[](0) == operator[](1) || operator[](1) == -1)
    {
        operator[](1) = -1;
        n--;
    }
    else if (operator[](1) == operator[](2) || operator[](2) == -1)
    {
        operator[](2) = -1;
        n--;
    }
    if (operator[](0) == operator[](2))
    {
        operator[](2) = -1;
        n--;
    }

    return n;
}


inline void Foam::triFace::flip()
{
    Swap(operator[](1), operator[](2));
}


inline Foam::pointField Foam::triFace::points(const UList<point>& points) const
{
    pointField p(3);

    p[0] = points[operator[](0)];
    p[1] = points[operator[](1)];
    p[2] = points[operator[](2)];

    return p;
}


inline Foam::face Foam::triFace::triFaceFace() const
{
    return Foam::face(*this);
}


inline Foam::triPointRef Foam::triFace::tri(const UList<point>& points) const
{
    return triPointRef
    (
        points[operator[](0)],
        points[operator[](1)],
        points[operator[](2)]
    );
}


inline Foam::point Foam::triFace::centre(const UList<point>& points) const
{
    return (1.0/3.0)*
    (
        points[operator[](0)]
      + points[operator[](1)]
      + points[operator[](2)]
    );
}


inline Foam::scalar Foam::triFace::mag(const UList<point>& points) const
{
    return ::Foam::mag(normal(points));
}


// could also delegate to triPointRef(...).normal()
inline Foam::vector Foam::triFace::normal(const UList<point>& points) const
{
    return 0.5*
    (
        (points[operator[](1)] - points[operator[](0)])
       ^(points[operator[](2)] - points[operator[](0)])
    );
}


inline Foam::label Foam::triFace::nTriangles() const
{
    return 1;
}


inline Foam::triFace Foam::triFace::reverseFace() const
{
    // The starting points of the original and reverse face are identical.
    return triFace(operator[](0), operator[](2), operator[](1));
}


inline bool Foam::triFace::found(const label globalIndex) const
{
    return which(globalIndex) != -1;
}


inline Foam::label Foam::triFace::which(const label globalIndex) const
{
    if (operator[](0) == globalIndex) return 0;
    if (operator[](1) == globalIndex) return 1;
    if (operator[](2) == globalIndex) return 2;
    return -1;
}


inline Foam::scalar Foam::triFace::sweptVol
(
    const UList<point>& opts,
    const UList<point>& npts
) const
{
    return (1.0/6.0)*
    (
        (
            (npts[operator[](0)] - opts[operator[](0)])
          & (
                (opts[operator[](1)] - opts[operator[](0)])
              ^ (opts[operator[](2)] - opts[operator[](0)])
            )
        )
      + (
            (npts[operator[](1)] - opts[operator[](1)])
          & (
                (opts[operator[](2)] - opts[operator[](1)])
              ^ (npts[operator[](0)] - opts[operator[](1)])
            )
        )
      + (
            (opts[operator[](2)] - npts[operator[](2)])
          & (
                (npts[operator[](1)] - npts[operator[](2)])
              ^ (npts[operator[](0)] - npts[operator[](2)])
            )
        )
    );
}


Foam::tensor Foam::triFace::inertia
(
    const UList<point>& points,
    const point& refPt,
    scalar density
) const
{
    // a triangle, do a direct calculation
    return this->tri(points).inertia(refPt, density);
}


inline Foam::pointHit Foam::triFace::ray
(
    const point& p,
    const vector& q,
    const UList<point>& points,
    const intersection::algorithm alg,
    const intersection::direction dir
) const
{
    return this->tri(points).ray(p, q, alg, dir);
}



inline Foam::pointHit Foam::triFace::intersection
(
    const point& p,
    const vector& q,
    const UList<point>& points,
    const intersection::algorithm alg,
    const scalar tol
) const
{
    return this->tri(points).intersection(p, q, alg, tol);
}


inline Foam::pointHit Foam::triFace::intersection
(
    const point& p,
    const vector& q,
    const point& ctr,
    const UList<point>& points,
    const intersection::algorithm alg,
    const scalar tol
) const
{
    return intersection(p, q, points, alg, tol);
}


inline Foam::pointHit Foam::triFace::nearestPoint
(
    const point& p,
    const UList<point>& points
) const
{
    return this->tri(points).nearestPoint(p);
}


inline Foam::pointHit Foam::triFace::nearestPointClassify
(
    const point& p,
    const UList<point>& points,
    label& nearType,
    label& nearLabel
) const
{
    return this->tri(points).nearestPointClassify(p, nearType, nearLabel);
}


inline Foam::label Foam::triFace::nEdges() const
{
    return 3;
}


inline Foam::edgeList Foam::triFace::edges() const
{
    edgeList e(3);

    e[0].first()  = operator[](0);
    e[0].second() = operator[](1);

    e[1].first()  = operator[](1);
    e[1].second() = operator[](2);

    e[2].first()  = operator[](2);
    e[2].second() = operator[](0);

    return e;
}


inline Foam::edge Foam::triFace::faceEdge(const label n) const
{
    return edge(operator[](n), operator[](fcIndex(n)));
}


// return
//  - +1: forward (counter-clockwise) on the face
//  - -1: reverse (clockwise) on the face
//  -  0: edge not found on the face
inline int Foam::triFace::edgeDirection(const edge& e) const
{
    if
    (
        (operator[](0) == e.first() && operator[](1) == e.second())
     || (operator[](1) == e.first() && operator[](2) == e.second())
     || (operator[](2) == e.first() && operator[](0) == e.second())
    )
    {
        return 1;
    }
    else if
    (
        (operator[](0) == e.second() && operator[](1) == e.first())
     || (operator[](1) == e.second() && operator[](2) == e.first())
     || (operator[](2) == e.second() && operator[](0) == e.first())
    )
    {
        return -1;
    }

    return 0;
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

inline bool Foam::operator==(const triFace& a, const triFace& b)
{
    return triFace::compare(a,b) != 0;
}


inline bool Foam::operator!=(const triFace& a, const triFace& b)
{
    return triFace::compare(a,b) == 0;
}


// ************************************************************************* //
