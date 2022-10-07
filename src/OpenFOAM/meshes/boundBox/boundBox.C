/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

#include "boundBox.H"
#include "PstreamReduceOps.H"
#include "plane.H"
#include "triangle.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::boundBox Foam::boundBox::greatBox
(
    point::uniform(-ROOTVGREAT),
    point::uniform(ROOTVGREAT)
);

const Foam::boundBox Foam::boundBox::invertedBox
(
    point::uniform(ROOTVGREAT),
    point::uniform(-ROOTVGREAT)
);

const Foam::faceList Foam::boundBox::faces
({
    // Point and face order as per hex cellmodel
    face({0, 4, 7, 3}),  // 0: x-min, left
    face({1, 2, 6, 5}),  // 1: x-max, right
    face({0, 1, 5, 4}),  // 2: y-min, bottom
    face({3, 7, 6, 2}),  // 3: y-max, top
    face({0, 3, 2, 1}),  // 4: z-min, back
    face({4, 5, 6, 7})   // 5: z-max, front
});

const Foam::FixedList<Foam::vector, 6> Foam::boundBox::faceNormals
({
    vector(-1,  0,  0), // 0: x-min, left
    vector( 1,  0,  0), // 1: x-max, right
    vector( 0, -1,  0), // 2: y-min, bottom
    vector( 0,  1,  0), // 3: y-max, top
    vector( 0,  0, -1), // 4: z-min, back
    vector( 0,  0,  1)  // 5: z-max, front
});


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::boundBox::boundBox(const boundBox& bb, const bool doReduce)
:
    boundBox(bb)
{
    if (doReduce)
    {
        reduce();
    }
}


Foam::boundBox::boundBox(const UList<point>& points, bool doReduce)
:
    boundBox()
{
    add(points);

    if (doReduce)
    {
        reduce();
    }
}


Foam::boundBox::boundBox(const tmp<pointField>& tpoints, bool doReduce)
:
    boundBox()
{
    add(tpoints);

    if (doReduce)
    {
        reduce();
    }
}


Foam::boundBox::boundBox
(
    const UList<point>& points,
    const labelUList& indices,
    bool doReduce
)
:
    boundBox()
{
    add(points, indices);

    if (doReduce)
    {
        reduce();
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::boundBox::points() const
{
    auto tpts = tmp<pointField>::New(8);
    auto& pts = tpts.ref();

    pts[0] = hexCorner<0>();
    pts[1] = hexCorner<1>();
    pts[2] = hexCorner<2>();
    pts[3] = hexCorner<3>();
    pts[4] = hexCorner<4>();
    pts[5] = hexCorner<5>();
    pts[6] = hexCorner<6>();
    pts[7] = hexCorner<7>();

    return tpts;
}


Foam::tmp<Foam::pointField> Foam::boundBox::faceCentres() const
{
    auto tpts = tmp<pointField>::New(6);
    auto& pts = tpts.ref();

    forAll(pts, facei)
    {
        pts[facei] = faceCentre(facei);
    }

    return tpts;
}


Foam::point Foam::boundBox::faceCentre(const direction facei) const
{
    point pt = boundBox::centre();

    switch (facei)
    {
        case 0: pt.x() = min().x(); break;  // 0: x-min, left
        case 1: pt.x() = max().x(); break;  // 1: x-max, right
        case 2: pt.y() = min().y(); break;  // 2: y-min, bottom
        case 3: pt.y() = max().y(); break;  // 3: y-max, top
        case 4: pt.z() = min().z(); break;  // 4: z-min, back
        case 5: pt.z() = max().z(); break;  // 5: z-max, front
        default:
        {
            FatalErrorInFunction
                << "Face:" << int(facei) << " should be [0..5]"
                << abort(FatalError);
        }
    }

    return pt;
}


void Foam::boundBox::reduce()
{
    Foam::reduce(min_, minOp<point>());
    Foam::reduce(max_, maxOp<point>());
}


Foam::boundBox Foam::boundBox::returnReduce(const boundBox& bb)
{
    boundBox work(bb);
    work.reduce();
    return work;
}


bool Foam::boundBox::intersect(const boundBox& bb)
{
    min_ = ::Foam::max(min_, bb.min_);
    max_ = ::Foam::min(max_, bb.max_);

    return valid();
}


bool Foam::boundBox::intersects(const plane& pln) const
{
    // Require a full 3D box
    if (nDim() != 3)
    {
        return false;
    }

    // Check as below(1) or above(2) - stop when it cuts both
    int side = 0;

    #undef  doLocalCode
    #define doLocalCode(Idx)                                                  \
    {                                                                         \
        side |= (pln.whichSide(hexCorner<Idx>()) == plane::BACK ? 1 : 2);     \
        if (side == 3) return true;  /* Both below and above: done */         \
    }

    // Each box corner
    doLocalCode(0);
    doLocalCode(1);
    doLocalCode(2);
    doLocalCode(3);
    doLocalCode(4);
    doLocalCode(5);
    doLocalCode(6);
    doLocalCode(7);

    #undef doLocalCode

    return false;
}


bool Foam::boundBox::contains(const UList<point>& points) const
{
    if (points.empty())
    {
        return true;
    }

    for (const point& p : points)
    {
        if (!contains(p))
        {
            return false;
        }
    }

    return true;
}


bool Foam::boundBox::containsAny(const UList<point>& points) const
{
    if (points.empty())
    {
        return true;
    }

    for (const point& p : points)
    {
        if (contains(p))
        {
            return true;
        }
    }

    return false;
}


Foam::point Foam::boundBox::nearest(const point& p) const
{
    // Clip the point to the range of the bounding box
    return point
    (
        Foam::min(Foam::max(p.x(), min_.x()), max_.x()),
        Foam::min(Foam::max(p.y(), min_.y()), max_.y()),
        Foam::min(Foam::max(p.z(), min_.z()), max_.z())
    );
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const boundBox& bb)
{
    if (os.format() == IOstreamOption::ASCII)
    {
        os << bb.min_ << token::SPACE << bb.max_;
    }
    else
    {
        os.write
        (
            reinterpret_cast<const char*>(&bb.min_),
            sizeof(boundBox)
        );
    }

    os.check(FUNCTION_NAME);
    return os;
}


Foam::Istream& Foam::operator>>(Istream& is, boundBox& bb)
{
    if (is.format() == IOstreamOption::ASCII)
    {
        is >> bb.min_ >> bb.max_;
    }
    else
    {
        Detail::readContiguous<boundBox>
        (
            is,
            reinterpret_cast<char*>(&bb.min_),
            sizeof(boundBox)
        );
    }

    is.check(FUNCTION_NAME);
    return is;
}


// ************************************************************************* //
