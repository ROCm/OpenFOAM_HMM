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

#include "scalarRange.H"
#include "string.H"
#include "error.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

bool Foam::scalarRange::parse(const std::string& str, scalarRange& range)
{
    range.clear();

    const auto colon = str.find(':');

    if (colon == std::string::npos)
    {
        // No colon

        if (str == "none")
        {
            // "none" is an empty (inverse) range
            return true;
        }

        // "VALUE"
        scalar val;
        if (readScalar(str, val))
        {
            range = scalarRange(val);
        }
    }
    else if (str.find(':', colon+1) != std::string::npos)
    {
        // A second colon is a syntax error
        return false;
    }
    else if (colon == 0)
    {
        // ":MAX"
        scalar val;
        if (readScalar(str.substr(1), val))
        {
            range = scalarRange::le(val);
        }
    }
    else if (colon == str.size()-1)
    {
        // "MIN:"
        scalar val;
        if (readScalar(str.substr(0, colon), val))
        {
            range = scalarRange::ge(val);
        }
    }
    else
    {
        // "MIN:MAX"
        scalar minVal, maxVal;
        if
        (
            readScalar(str.substr(0, colon), minVal)
         && readScalar(str.substr(colon+1), maxVal)
        )
        {
            range = scalarRange(minVal, maxVal);
        }
    }

    return range.valid();
}


Foam::scalarRange Foam::scalarRange::parse(const std::string& str)
{
    scalarRange range;

    if (!parse(str, range))
    {
        Info<< "Bad scalar-range while parsing: " << str << endl;
    }

    return range;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const scalarRange& range)
{
    switch (range.type_)
    {
        case scalarRange::EQ:
            os << range.min_;
            break;

        case scalarRange::GE:
        case scalarRange::GT:
            os << range.min_ << ":Inf";
            break;

        case scalarRange::LE:
        case scalarRange::LT:
            os << "-Inf:" << range.max_;
            break;

        case scalarRange::GE_LE:
            os << range.min_ << ':' << range.max_;
            break;

        default:
            os << "none";
            break;
    }

    return os;
}


// ************************************************************************* //
