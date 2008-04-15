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

\*---------------------------------------------------------------------------*/

#include "IOstreams.H"
#include "face.H"
#include "triPointRef.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct null
inline triFace::triFace()
{}


//- Construct from components
inline triFace::triFace
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


inline triFace::triFace(Istream& is)
:
    FixedList<label, 3>(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return the points associated with this face
inline pointField triFace::points(const pointField& points) const
{
    pointField p(3);

    p[0] = points[operator[](0)];
    p[1] = points[operator[](1)];
    p[2] = points[operator[](2)];

    return p;
}


// Return triagle as a face
inline face triFace::triFaceFace() const
{
    face f(3);

    f[0] = operator[](0);
    f[1] = operator[](1);
    f[2] = operator[](2);

    return f;
}


inline label triFace::nEdges() const
{
    return 3;
}


inline edgeList triFace::edges() const
{
    edgeList e(3);

    e[0].start() = operator[](0);
    e[0].end() = operator[](1);

    e[1].start() = operator[](1);
    e[1].end() = operator[](2);

    e[2].start() = operator[](2);
    e[2].end() = operator[](0);

    return e;
}


inline point triFace::centre(const pointField& points) const
{
    return (1.0/3.0)*
    (
        points[operator[](0)]
      + points[operator[](1)]
      + points[operator[](2)]
    );
}


inline scalar triFace::mag(const pointField& points) const
{
    return ::Foam::mag(normal(points));
}


inline vector triFace::normal(const pointField& points) const
{
    return 0.5*
    (
        (points[operator[](1)] - points[operator[](0)])
       ^(points[operator[](2)] - points[operator[](0)])
    );
}


inline scalar triFace::sweptVol
(
    const pointField& opts,
    const pointField& npts
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


inline pointHit triFace::ray
(
    const point& p,
    const vector& q,
    const pointField& points,
    const intersection::algorithm alg,
    const intersection::direction dir
) const
{
    return triPointRef
    (
        points[operator[](0)],
        points[operator[](1)],
        points[operator[](2)]
    ).ray(p, q, alg, dir);
}


inline triPointRef triFace::tri(const pointField& points) const
{
    return triPointRef
    (
        points[operator[](0)],
        points[operator[](1)],
        points[operator[](2)]
    );
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline bool operator==(const triFace& tf1, const triFace& tf2)
{
    return
    (
        (tf1[0] == tf2[0] && tf1[1] == tf2[1] && tf1[2] == tf2[2])
     || (tf1[0] == tf2[1] && tf1[1] == tf2[2] && tf1[2] == tf2[0])
     || (tf1[0] == tf2[2] && tf1[1] == tf2[0] && tf1[2] == tf2[1])
     || (tf1[0] == tf2[2] && tf1[1] == tf2[1] && tf1[2] == tf2[0])
     || (tf1[0] == tf2[1] && tf1[1] == tf2[0] && tf1[2] == tf2[2])
     || (tf1[0] == tf2[0] && tf1[1] == tf2[2] && tf1[2] == tf2[1])
    );
}


inline bool operator!=(const triFace& tf1, const triFace& tf2)
{
    return !(tf1 == tf2);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
