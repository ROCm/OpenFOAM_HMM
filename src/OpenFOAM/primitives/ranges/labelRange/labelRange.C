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

int Foam::labelRange::debug(::Foam::debug::debugSwitch("labelRange", 0));


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
        size_ += start_;
        start_ = 0;
    }

    if (size_ < 0)
    {
        size_ = 0;
    }
}


bool Foam::labelRange::overlaps(const labelRange& range, bool touches) const
{
    const label final = touches ? 1 : 0;

    return
    (
        this->size()
     && range.size()
     &&
        (
            (
                range.first() >= this->first()
             && range.first() <= this->last() + final
            )
         ||
            (
                this->first() >= range.first()
             && this->first() <= range.last() + final
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
    else if (!range.size_)
    {
        return range;
    }

    const label lower = Foam::min(this->first(), range.first());
    const label upper = Foam::max(this->last(),  range.last());
    const label sz = upper - lower + 1;

    return labelRange(lower, sz);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::labelRange::operator+=(const labelRange& rhs)
{
    if (!size_)
    {
        // trivial case
        operator=(rhs);
    }
    else if (rhs.size_)
    {
        const label lower = Foam::min(this->first(), rhs.first());
        const label upper = Foam::max(this->last(),  rhs.last());

        start_ = lower;
        size_  = upper - lower + 1;
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, labelRange& range)
{
    is.readBegin("labelRange");
    is  >> range.start_ >> range.size_;
    is.readEnd("labelRange");

    is.check("operator>>(Istream&, labelRange&)");

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
        << range.start_ << token::SPACE << range.size_
        << token::END_LIST;

    os.check("operator<<(Ostream&, const labelRange&)");
    return os;
}


// ************************************************************************* //
