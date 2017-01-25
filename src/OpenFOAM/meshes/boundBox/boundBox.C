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

\*---------------------------------------------------------------------------*/

#include "boundBox.H"
#include "PstreamReduceOps.H"
#include "tmp.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// (min,max) = (-VGREAT,+VGREAT)
const Foam::boundBox Foam::boundBox::greatBox(point::min, point::max);

// (min,max) = (+VGREAT,-VGREAT)
const Foam::boundBox Foam::boundBox::invertedBox(point::max, point::min);

const Foam::faceList Foam::boundBox::faces
({
    // Point and face order as per hex cellmodel
    face{0, 4, 7, 3}, // x-min
    face{1, 2, 6, 5}, // x-max
    face{0, 1, 5, 4}, // y-min
    face{3, 7, 6, 2}, // y-max
    face{0, 3, 2, 1}, // z-min
    face{4, 5, 6, 7}  // z-max
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::boundBox::calculate(const UList<point>& points, const bool doReduce)
{
    if (points.empty())
    {
        min_ = Zero;
        max_ = Zero;

        if (doReduce && Pstream::parRun())
        {
            // Use values that get overwritten by reduce minOp, maxOp below
            min_ = point(VGREAT, VGREAT, VGREAT);
            max_ = point(-VGREAT, -VGREAT, -VGREAT);
        }
    }
    else
    {
        min_ = points[0];
        max_ = points[0];

        for (label i = 1; i < points.size(); i++)
        {
            min_ = ::Foam::min(min_, points[i]);
            max_ = ::Foam::max(max_, points[i]);
        }
    }

    // Reduce parallel information
    if (doReduce)
    {
        reduce();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::boundBox::boundBox(const UList<point>& points, const bool doReduce)
:
    min_(Zero),
    max_(Zero)
{
    calculate(points, doReduce);
}


Foam::boundBox::boundBox(const tmp<pointField>& points, const bool doReduce)
:
    min_(Zero),
    max_(Zero)
{
    calculate(points(), doReduce);
    points.clear();
}


Foam::boundBox::boundBox
(
    const UList<point>& points,
    const labelUList& indices,
    const bool doReduce
)
:
    min_(Zero),
    max_(Zero)
{
    if (points.empty() || indices.empty())
    {
        if (doReduce && Pstream::parRun())
        {
            // Use values that get overwritten by reduce minOp, maxOp below
            min_ = point(VGREAT, VGREAT, VGREAT);
            max_ = point(-VGREAT, -VGREAT, -VGREAT);
        }
    }
    else
    {
        min_ = points[indices[0]];
        max_ = points[indices[0]];

        for (label i=1; i < indices.size(); ++i)
        {
            min_ = ::Foam::min(min_, points[indices[i]]);
            max_ = ::Foam::max(max_, points[indices[i]]);
        }
    }

    // Reduce parallel information
    if (doReduce)
    {
        reduce();
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::boundBox::points() const
{
    tmp<pointField> tPts = tmp<pointField>(new pointField(8));
    pointField& pt = tPts.ref();

    pt[0] = min_;                                   // min-x, min-y, min-z
    pt[1] = point(max_.x(), min_.y(), min_.z());    // max-x, min-y, min-z
    pt[2] = point(max_.x(), max_.y(), min_.z());    // max-x, max-y, min-z
    pt[3] = point(min_.x(), max_.y(), min_.z());    // min-x, max-y, min-z
    pt[4] = point(min_.x(), min_.y(), max_.z());    // min-x, min-y, max-z
    pt[5] = point(max_.x(), min_.y(), max_.z());    // max-x, min-y, max-z
    pt[6] = max_;                                   // max-x, max-y, max-z
    pt[7] = point(min_.x(), max_.y(), max_.z());    // min-x, max-y, max-z

    return tPts;
}


void Foam::boundBox::inflate(const scalar s)
{
    vector ext = vector::one*s*mag();

    min_ -= ext;
    max_ += ext;
}


void Foam::boundBox::reduce()
{
    Foam::reduce(min_, minOp<point>());
    Foam::reduce(max_, maxOp<point>());
}


bool Foam::boundBox::contains(const UList<point>& points) const
{
    if (points.empty())
    {
        return true;
    }

    forAll(points, i)
    {
        if (!contains(points[i]))
        {
            return false;
        }
    }

    return true;
}


bool Foam::boundBox::contains
(
    const UList<point>& points,
    const labelUList& indices
) const
{
    if (points.empty() || indices.empty())
    {
        return true;
    }

    forAll(indices, i)
    {
        if (!contains(points[indices[i]]))
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

    forAll(points, i)
    {
        if (contains(points[i]))
        {
            return true;
        }
    }

    return false;
}


bool Foam::boundBox::containsAny
(
    const UList<point>& points,
    const labelUList& indices
) const
{
    if (points.empty() || indices.empty())
    {
        return true;
    }

    forAll(indices, i)
    {
        if (contains(points[indices[i]]))
        {
            return true;
        }
    }

    return false;
}


Foam::point Foam::boundBox::nearest(const point& pt) const
{
    // Clip the point to the range of the bounding box
    const scalar surfPtx = Foam::max(Foam::min(pt.x(), max_.x()), min_.x());
    const scalar surfPty = Foam::max(Foam::min(pt.y(), max_.y()), min_.y());
    const scalar surfPtz = Foam::max(Foam::min(pt.z(), max_.z()), min_.z());

    return point(surfPtx, surfPty, surfPtz);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const boundBox& bb)
{
    if (os.format() == IOstream::ASCII)
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

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const boundBox&)");
    return os;
}


Foam::Istream& Foam::operator>>(Istream& is, boundBox& bb)
{
    if (is.format() == IOstream::ASCII)
    {
        is >> bb.min_ >> bb.max_;
    }
    else
    {
        is.read
        (
            reinterpret_cast<char*>(&bb.min_),
            sizeof(boundBox)
        );
    }

    // Check state of Istream
    is.check("Istream& operator>>(Istream&, boundBox&)");
    return is;
}


// ************************************************************************* //
