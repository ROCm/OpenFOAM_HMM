/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
     \\/     M anipulation  |
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::scalarRange::scalarRange
(
    const scalarRange::testType type,
    const scalar minVal,
    const scalar maxVal
)
:
    type_(type),
    min_(minVal),
    max_(maxVal)
{}


inline Foam::scalarRange::scalarRange()
:
    scalarRange(NONE, Foam::VGREAT, -Foam::VGREAT)
{}


inline Foam::scalarRange::scalarRange(const scalar value)
:
    scalarRange(EQ, value, value)
{}


inline Foam::scalarRange::scalarRange(const scalar minVal, const scalar maxVal)
:
    scalarRange(GE_LE, minVal, maxVal)
{
    if (minVal > maxVal)
    {
        clear();  // Inverted - explicitly mark as such
    }
    else if (minVal == maxVal)
    {
        type_ = EQ;
    }
}


inline Foam::scalarRange Foam::scalarRange::ge(const scalar minVal)
{
    return scalarRange(GE, minVal, Foam::VGREAT);
}


inline Foam::scalarRange Foam::scalarRange::gt(const scalar minVal)
{
    return scalarRange(GT, minVal, Foam::VGREAT);
}


inline Foam::scalarRange Foam::scalarRange::le(const scalar maxVal)
{
    return scalarRange(LE, -Foam::VGREAT, maxVal);
}

inline Foam::scalarRange Foam::scalarRange::lt(const scalar maxVal)
{
    return scalarRange(LT, -Foam::VGREAT, maxVal);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::scalarRange::clear()
{
    type_ = scalarRange::NONE;
    min_ = Foam::VGREAT;
    max_ = -Foam::VGREAT;
}


inline bool Foam::scalarRange::empty() const
{
    return type_ == NONE;
}


inline bool Foam::scalarRange::valid() const
{
    return type_;
}


inline bool Foam::scalarRange::single() const
{
    return type_ == EQ;
}


inline Foam::scalar Foam::scalarRange::min() const
{
    return min_;
}


inline Foam::scalar Foam::scalarRange::max() const
{
    return max_;
}


inline Foam::scalar Foam::scalarRange::value() const
{
    switch (type_)
    {
        case EQ:
            return min_;    // For equals, min and max are identical

        case GE:
        case GT:
            return min_;

        case LE:
        case LT:
            return max_;

        case GE_LE:
            return 0.5*(min_ + max_);

        default:
            return 0;
    }
}


inline bool Foam::scalarRange::match(const scalar& value) const
{
    switch (type_)
    {
        case EQ:    return equal(value, min_);
        case GE:    return (value >= min_);
        case GT:    return (value > min_);
        case LE:    return (value <= max_);
        case LT:    return (value < max_);
        case GE_LE: return (value >= min_ && value <= max_);
        default:    return false;
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::scalarRange::operator()(const scalar& value) const
{
    return match(value);
}


inline bool Foam::scalarRange::operator==(const scalarRange& rhs) const
{
    return (type_ == rhs.type_ && min_ == rhs.min_ && max_ == rhs.max_);
}


inline bool Foam::scalarRange::operator!=(const scalarRange& rhs) const
{
    return !(operator==(rhs));
}


// ************************************************************************* //
