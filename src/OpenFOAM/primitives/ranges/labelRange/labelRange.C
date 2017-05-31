/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "labelRange.H"
#include "token.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
int labelRange::debug(debug::debugSwitch("labelRange", 0));
}

const Foam::labelRange Foam::labelRange::null;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::labelRange::labelRange(Istream& is)
:
    start_(0),
    size_(0)
{
    is  >> *this;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::labelRange::adjust()
{
    if (start_ < 0)
    {
        if (size_ <= 0)
        {
            size_ = 0;
        }
        else
        {
            size_ += start_;
        }
        start_ = 0;
    }
    else if (size_ < 0)
    {
        size_ = 0;
    }
}


bool Foam::labelRange::overlaps(const labelRange& range, bool touches) const
{
    const label extra = touches ? 1 : 0;

    return
    (
        this->size() && range.size()
     &&
        (
            (
                range.first() >= this->first()
             && range.first() <= this->last() + extra
            )
         ||
            (
                this->first() >= range.first()
             && this->first() <= range.last() + extra
            )
        )
    );
}


Foam::labelRange Foam::labelRange::join(const labelRange& range) const
{
    // trivial cases first
    if (!size_)
    {
        return *this;
    }
    else if (!range.size())
    {
        return range;
    }

    const label lower = Foam::min(this->first(), range.first());
    const label upper = Foam::max(this->last(),  range.last());
    const label total = upper+1 - lower;
    // last = start+size-1
    // size = last+1-start

    return labelRange(lower, total);
}


Foam::labelRange Foam::labelRange::subset(const labelRange& range) const
{
    const label lower = Foam::max(this->first(), range.first());
    const label upper = Foam::min(this->last(),  range.last());
    const label total = upper+1 - lower;
    // last = start+size-1
    // size = last+1-start

    if (total > 0)
    {
        return labelRange(lower, total);
    }
    else
    {
        return labelRange();
    }
}


Foam::labelRange Foam::labelRange::subset
(
    const label start,
    const label size
) const
{
    const label lower = Foam::max(this->start(), start);
    const label upper = Foam::min(this->last(),  start+Foam::max(0,size-1));
    const label total = upper+1 - lower;
    // last = start+size-1
    // size = last+1-start

    if (total > 0)
    {
        return labelRange(lower, total);
    }
    else
    {
        return labelRange();
    }
}


Foam::labelRange Foam::labelRange::subset0(const label size) const
{
    const label lower = Foam::max(this->start(), 0);
    const label upper = Foam::min(this->last(),  Foam::max(0,size-1));
    const label total = upper+1 - lower;
    // last = start+size-1
    // size = last+1-start

    if (total > 0)
    {
        return labelRange(lower, total);
    }
    else
    {
        return labelRange();
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, labelRange& range)
{
    is.readBegin("labelRange");
    is  >> range.start_ >> range.size_;
    is.readEnd("labelRange");

    is.check(FUNCTION_NAME);

    // Disallow invalid sizes
    if (range.size_ < 0)
    {
        range.size_ = 0;
    }

    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const labelRange& range)
{
    // Write ASCII only for now
    os  << token::BEGIN_LIST
        << range.start() << token::SPACE << range.size()
        << token::END_LIST;

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
