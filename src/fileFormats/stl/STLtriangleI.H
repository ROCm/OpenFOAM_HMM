/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2016 OpenCFD Ltd.
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

#include "triPointRef.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::STLtriangle::STLtriangle()
{}


inline Foam::STLtriangle::STLtriangle
(
    const STLpoint& normal,
    const STLpoint& a,
    const STLpoint& b,
    const STLpoint& c,
    unsigned short attrib
)
:
    normal_(normal),
    a_(a),
    b_(b),
    c_(c),
    attrib_(attrib)
{}


inline Foam::STLtriangle::STLtriangle(std::istream& is)
{
    read(is);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::STLpoint& Foam::STLtriangle::normal() const
{
    return normal_;
}


inline const Foam::STLpoint& Foam::STLtriangle::a() const
{
    return a_;
}


inline const Foam::STLpoint& Foam::STLtriangle::b() const
{
    return b_;
}


inline const Foam::STLpoint& Foam::STLtriangle::c() const
{
    return c_;
}


inline unsigned short Foam::STLtriangle::attrib() const
{
    return attrib_;
}


inline void Foam::STLtriangle::read(std::istream& is)
{
    is.read(reinterpret_cast<char*>(&normal_), 4*sizeof(STLpoint));
    is.read(reinterpret_cast<char*>(&attrib_), sizeof(STLattrib));
}


inline void Foam::STLtriangle::write(std::ostream& os) const
{
    os.write(reinterpret_cast<const char*>(&normal_), 4*sizeof(STLpoint));
    os.write(reinterpret_cast<const char*>(&attrib_), sizeof(STLattrib));
}


inline Foam::Ostream& Foam::STLtriangle::print(Ostream& os) const
{
    os  << " facet normal "
        << normal_.x() << ' ' << normal_.y() << ' ' << normal_.z() << nl
        << "  outer loop" << nl
        << "   vertex " << a_.x() << ' ' << a_.y() << ' ' << a_.z() << nl
        << "   vertex " << b_.x() << ' ' << b_.y() << ' ' << b_.z() << nl
        << "   vertex " << c_.x() << ' ' << c_.y() << ' ' << c_.z() << nl
        << "  endloop" << nl
        << " endfacet" << nl;

    return os;
}


inline void Foam::STLtriangle::write
(
    Ostream& os,
    const vector& norm,
    const point& pt0,
    const point& pt1,
    const point& pt2
)
{
    os  << " facet normal "
        << norm.x() << ' ' << norm.y() << ' ' << norm.z() << nl
        << "  outer loop" << nl
        << "   vertex " << pt0.x() << ' ' << pt0.y() << ' ' << pt0.z() << nl
        << "   vertex " << pt1.x() << ' ' << pt1.y() << ' ' << pt1.z() << nl
        << "   vertex " << pt2.x() << ' ' << pt2.y() << ' ' << pt2.z() << nl
        << "  endloop" << nl
        << " endfacet" << nl;
}


inline void Foam::STLtriangle::write
(
    Ostream& os,
    const point& pt0,
    const point& pt1,
    const point& pt2
)
{
    // Calculate the normal ourselves
    const vector norm = triPointRef(pt0, pt1, pt2).unitNormal();

    write(os, norm, pt0, pt1, pt2);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

inline Foam::Ostream& Foam::operator<<(Ostream& os, const STLtriangle& tri)
{
    os  << tri.normal_ << token::SPACE
        << tri.a_ << token::SPACE
        << tri.b_ << token::SPACE
        << tri.c_ << token::SPACE
        << tri.attrib_;

    return os;
}


// ************************************************************************* //
