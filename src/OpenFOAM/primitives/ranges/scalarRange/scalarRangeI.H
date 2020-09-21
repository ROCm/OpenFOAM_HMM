/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2020 OpenCFD Ltd.
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

inline constexpr Foam::scalarRange::scalarRange
(
    const scalarRange::testType type,
    const scalar minVal,
    const scalar maxVal
) noexcept
:
    min_(minVal),
    max_(maxVal),
    type_(type)
{}


inline constexpr Foam::scalarRange::scalarRange() noexcept
:
    scalarRange(scalarRange::NONE, GREAT, -GREAT)
{}


inline constexpr Foam::scalarRange::scalarRange(const scalar val) noexcept
:
    scalarRange(scalarRange::EQ, val, val)
{}


inline Foam::scalarRange::scalarRange
(
    const scalar minVal,
    const scalar maxVal
) noexcept
:
    scalarRange(scalarRange::GE_LE, minVal, maxVal)
{
    if (minVal > maxVal)
    {
        clear();  // Inverted - explicitly mark as such
    }
    else if (equal(minVal, maxVal))
    {
        type_ = scalarRange::EQ;
    }
}


inline constexpr Foam::scalarRange
Foam::scalarRange::ge(const scalar minVal) noexcept
{
    return scalarRange(scalarRange::GE, minVal, GREAT);
}


inline constexpr Foam::scalarRange
Foam::scalarRange::gt(const scalar minVal) noexcept
{
    return scalarRange(scalarRange::GT, minVal, GREAT);
}


inline constexpr Foam::scalarRange
Foam::scalarRange::ge0() noexcept
{
    return scalarRange(scalarRange::GE, 0, GREAT);
}


inline constexpr Foam::scalarRange
Foam::scalarRange::gt0() noexcept
{
    return scalarRange(scalarRange::GT, 0, GREAT);
}


inline constexpr Foam::scalarRange
Foam::scalarRange::le(const scalar maxVal) noexcept
{
    return scalarRange(scalarRange::LE, -GREAT, maxVal);
}

inline constexpr Foam::scalarRange
Foam::scalarRange::lt(const scalar maxVal) noexcept
{
    return scalarRange(scalarRange::LT, -GREAT, maxVal);
}


inline constexpr Foam::scalarRange
Foam::scalarRange::zero_one() noexcept
{
    return scalarRange(scalarRange::GE_LE, 0, 1);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::scalarRange::clear() noexcept
{
    min_ = GREAT;
    max_ = -GREAT;
    type_ = scalarRange::NONE;
}


inline bool Foam::scalarRange::empty() const noexcept
{
    return type_ == scalarRange::NONE;
}


inline bool Foam::scalarRange::valid() const noexcept
{
    return bool(type_);
}


inline bool Foam::scalarRange::single() const noexcept
{
    return type_ == scalarRange::EQ;
}


inline Foam::scalar Foam::scalarRange::min() const noexcept
{
    return min_;
}


inline Foam::scalar Foam::scalarRange::max() const noexcept
{
    return max_;
}


inline Foam::scalar Foam::scalarRange::value() const
{
    switch (type_)
    {
        case scalarRange::EQ: // For equals, min and max are identical
        case scalarRange::GE:
        case scalarRange::GT:
            return min_;

        case scalarRange::LE:
        case scalarRange::LT:
            return max_;

        case scalarRange::GE_LE:
            // Multiply before adding to avoid possible overflow
            return (0.5 * min_) + (0.5 * max_);

        default:  // NONE or ALWAYS. Zero is reasonable dummy value.
            return 0;
    }
}


inline bool Foam::scalarRange::match(const scalar& val) const
{
    switch (type_)
    {
        case scalarRange::EQ:    return equal(val, min_);
        case scalarRange::GE:    return (val >= min_);
        case scalarRange::GT:    return (val > min_);
        case scalarRange::LE:    return (val <= max_);
        case scalarRange::LT:    return (val < max_);
        case scalarRange::GE_LE: return (val >= min_ && val <= max_);
        case scalarRange::ALWAYS: return true;
        default:    return false;
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::scalarRange::operator()(const scalar& val) const
{
    return match(val);
}


inline constexpr bool
Foam::scalarRange::operator==
(
    const scalarRange& rhs
) const noexcept
{
    return (type_ == rhs.type_ && min_ == rhs.min_ && max_ == rhs.max_);
}


inline constexpr bool
Foam::scalarRange::operator!=
(
    const scalarRange& rhs
) const noexcept
{
    return !(*this == rhs);
}


// ************************************************************************* //
