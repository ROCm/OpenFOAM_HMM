/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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

#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::objectMap::objectMap()
:
    index_(-1),
    masterObjects_()
{}


inline Foam::objectMap::objectMap(const label index, const UList<label>& master)
:
    index_(index),
    masterObjects_(master)
{}


inline Foam::objectMap::objectMap(Istream& is)
{
    is.readBegin("objectMap");

    is >> index_ >> static_cast<labelList&>(masterObjects_);

    is.readEnd("objectMap");
    is.check(FUNCTION_NAME);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label& Foam::objectMap::index()
{
    return index_;
}


inline Foam::label Foam::objectMap::index() const
{
    return index_;
}


inline Foam::labelList& Foam::objectMap::masterObjects()
{
    return masterObjects_;
}


inline const Foam::labelList& Foam::objectMap::masterObjects() const
{
    return masterObjects_;
}


// * * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * //

inline bool Foam::operator==(const objectMap& a, const objectMap& b)
{
    return
    (
        (a.index_ == b.index_) && (a.masterObjects_ == b.masterObjects_)
    );
}


inline bool Foam::operator!=(const objectMap& a, const objectMap& b)
{
    return !(a == b);
}


// * * * * * * * * * * * * * * * Ostream Operator *  * * * * * * * * * * * * //

inline Foam::Ostream& Foam::operator<<(Ostream& os, const objectMap& a)
{
    os  << token::BEGIN_LIST
        << a.index_ << token::SPACE
        << a.masterObjects_
        << token::END_LIST;

    os.check(FUNCTION_NAME);
    return os;
}


inline Foam::Istream& Foam::operator>>(Istream& is, objectMap& a)
{
    is.readBegin("objectMap");
    is >> a.index_ >> a.masterObjects_;
    is.readEnd("objectMap");

    is.check(FUNCTION_NAME);
    return is;
}


// ************************************************************************* //
