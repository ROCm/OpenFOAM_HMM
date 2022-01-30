/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2018 OpenCFD Ltd.
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

#include "coordSet.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::coordSet::coordFormat
>
Foam::coordSet::coordFormatNames
({
    { coordFormat::XYZ, "xyz" },
    { coordFormat::X, "x" },
    { coordFormat::Y, "y" },
    { coordFormat::Z, "z" },
    { coordFormat::DISTANCE, "distance" },
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::coordSet::checkDimensions() const
{
    if (size() != curveDist_.size())
    {
        FatalErrorInFunction
            << "Size of points and curve distance must be the same" << nl
            << "    points size : " << size()
            << "    curve size  : " << curveDist_.size()
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coordSet::coordSet
(
    const word& name,
    const coordFormat axisType
)
:
    pointField(),
    name_(name),
    axis_(axisType),
    curveDist_()
{}


Foam::coordSet::coordSet
(
    const word& name,
    const word& axis
)
:
    pointField(),
    name_(name),
    axis_(coordFormatNames[axis]),
    curveDist_()
{}


Foam::coordSet::coordSet
(
    const word& name,
    const word& axis,
    const List<point>& points,
    const scalarList& curveDist
)
:
    pointField(points),
    name_(name),
    axis_(coordFormatNames[axis]),
    curveDist_(curveDist)
{
    checkDimensions();
}


Foam::coordSet::coordSet
(
    const word& name,
    const word& axis,
    List<point>&& points,
    scalarList&& curveDist
)
:
    pointField(std::move(points)),
    name_(name),
    axis_(coordFormatNames[axis]),
    curveDist_(std::move(curveDist))
{
    checkDimensions();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::coordSet::hasVectorAxis() const
{
    return axis_ == coordFormat::XYZ;
}


Foam::scalar Foam::coordSet::scalarCoord(const label index) const
{
    const point& p = operator[](index);

    switch (axis_)
    {
        case coordFormat::X:
        {
            return p.x();
        }
        case coordFormat::Y:
        {
            return p.y();
        }
        case coordFormat::Z:
        {
            return p.z();
        }
        case coordFormat::DISTANCE:
        {
            // Note: If this has been constructed from the 'name' and 'axis'
            // constructor the curveDist list will not have been set

            if (curveDist_.empty())
            {
                FatalErrorInFunction
                    << "Axis type '" << coordFormatNames[axis_]
                    << "' requested but curve distance has not been set"
                    << abort(FatalError);
            }

            return curveDist_[index];
        }
        default:
        {
            FatalErrorInFunction
                << "Illegal axis specification '" << coordFormatNames[axis_]
                << "' for sampling line " << name_
                << exit(FatalError);

            return 0;
        }
    }
}


Foam::point Foam::coordSet::vectorCoord(const label index) const
{
    const point& p = operator[](index);

    return p;
}


Foam::Ostream& Foam::coordSet::write(Ostream& os) const
{
    os  << "name:" << name_ << " axis:" << coordFormatNames[axis_]
        << nl
        << nl << "\t(coord)"
        << endl;

    for (const point& pt : *this)
    {
        os  << '\t' << pt << endl;
    }

    return os;
}


// ************************************************************************* //
