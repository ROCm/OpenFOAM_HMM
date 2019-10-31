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

\*---------------------------------------------------------------------------*/

#include "boundBox.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::boundBox::boundBox()
:
    min_(invertedBox.min()),
    max_(invertedBox.max())
{}


inline Foam::boundBox::boundBox(const point& pt)
:
    min_(pt),
    max_(pt)
{}


inline Foam::boundBox::boundBox(const point& min, const point& max)
:
    min_(min),
    max_(max)
{}


inline Foam::boundBox::boundBox(Istream& is)
{
    operator>>(is, *this);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::boundBox::empty() const
{
    // Check each component for max < min
    for (direction dir = 0; dir < vector::nComponents; ++dir)
    {
        if (max_[dir] < min_[dir])
        {
            return true;
        }
    }
    return false;
}


inline bool Foam::boundBox::valid() const
{
    // Check each component for max < min
    for (direction dir = 0; dir < vector::nComponents; ++dir)
    {
        if (max_[dir] < min_[dir])
        {
            return false;
        }
    }

    return true;
}


inline const Foam::point& Foam::boundBox::min() const
{
    return min_;
}


inline const Foam::point& Foam::boundBox::max() const
{
    return max_;
}


inline Foam::point& Foam::boundBox::min()
{
    return min_;
}


inline Foam::point& Foam::boundBox::max()
{
    return max_;
}


inline Foam::point Foam::boundBox::centre() const
{
    return 0.5 * (min_ + max_);
}


inline Foam::point Foam::boundBox::midpoint() const
{
    return this->centre();
}


inline Foam::vector Foam::boundBox::span() const
{
    return (max_ - min_);
}


inline Foam::scalar Foam::boundBox::mag() const
{
    return ::Foam::mag(max_ - min_);
}


inline Foam::scalar Foam::boundBox::volume() const
{
    return cmptProduct(span());
}


inline Foam::scalar Foam::boundBox::minDim() const
{
    return cmptMin(span());
}


inline Foam::scalar Foam::boundBox::maxDim() const
{
    return cmptMax(span());
}


inline Foam::scalar Foam::boundBox::avgDim() const
{
    return cmptAv(span());
}


inline Foam::label Foam::boundBox::nDim() const
{
    label ngood = 0;

    for (direction dir = 0; dir < vector::nComponents; ++dir)
    {
        const scalar diff = (max_[dir] - min_[dir]);
        if (diff < 0)
        {
            return -1;
        }
        else if (diff > 0)
        {
            ++ngood;
        }
    }

    return ngood;
}


inline void Foam::boundBox::clear()
{
    min_ = invertedBox.min();
    max_ = invertedBox.max();
}


inline void Foam::boundBox::add(const boundBox& bb)
{
    min_ = ::Foam::min(min_, bb.min_);
    max_ = ::Foam::max(max_, bb.max_);
}


inline void Foam::boundBox::add(const point& pt)
{
    min_ = ::Foam::min(min_, pt);
    max_ = ::Foam::max(max_, pt);
}


inline void Foam::boundBox::add(const UList<point>& points)
{
    for (const point& p : points)
    {
        add(p);
    }
}


inline void Foam::boundBox::add(const tmp<pointField>& tpoints)
{
    add(tpoints());
    tpoints.clear();
}


inline bool Foam::boundBox::overlaps(const boundBox& bb) const
{
    return
    (
        bb.max_.x() >= min_.x() && bb.min_.x() <= max_.x()
     && bb.max_.y() >= min_.y() && bb.min_.y() <= max_.y()
     && bb.max_.z() >= min_.z() && bb.min_.z() <= max_.z()
    );
}


inline bool Foam::boundBox::overlaps
(
    const point& centre,
    const scalar radiusSqr
) const
{
    // Find out where centre is in relation to bb.
    // Find nearest point on bb.
    scalar distSqr = 0;

    for (direction dir = 0; dir < vector::nComponents; ++dir)
    {
        const scalar d0 = min_[dir] - centre[dir];
        const scalar d1 = max_[dir] - centre[dir];

        if ((d0 > 0) != (d1 > 0))
        {
            // centre inside both extrema. This component does not add any
            // distance.
        }
        else if (Foam::mag(d0) < Foam::mag(d1))
        {
            distSqr += d0*d0;
        }
        else
        {
            distSqr += d1*d1;
        }

        if (distSqr > radiusSqr)
        {
            return false;
        }
    }

    return true;
}


inline bool Foam::boundBox::contains(const point& pt) const
{
    return
    (
        min_.x() <= pt.x() && pt.x() <= max_.x()
     && min_.y() <= pt.y() && pt.y() <= max_.y()
     && min_.z() <= pt.z() && pt.z() <= max_.z()
    );
}


inline bool Foam::boundBox::contains(const boundBox& bb) const
{
    return contains(bb.min()) && contains(bb.max());
}


inline bool Foam::boundBox::containsInside(const point& pt) const
{
    return
    (
        min_.x() < pt.x() && pt.x() < max_.x()
     && min_.y() < pt.y() && pt.y() < max_.y()
     && min_.z() < pt.z() && pt.z() < max_.z()
    );
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void Foam::boundBox::operator+=(const boundBox& bb)
{
    add(bb);
}


// * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * * //

inline bool Foam::operator==(const boundBox& a, const boundBox& b)
{
    return (a.min() == b.min()) && (a.max() == b.max());
}


inline bool Foam::operator!=(const boundBox& a, const boundBox& b)
{
    return !(a == b);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
