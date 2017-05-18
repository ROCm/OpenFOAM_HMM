/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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

#include "error.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const char* Foam::ensightCells::key(const enum elemType what)
{
    return elemNames[what];
}


inline Foam::label Foam::ensightCells::index() const
{
    return index_;
}


inline Foam::label& Foam::ensightCells::index()
{
    return index_;
}


inline Foam::label Foam::ensightCells::size() const
{
    return address_.size();
}


inline const Foam::FixedList<Foam::label,5>& Foam::ensightCells::totals() const
{
    return sizes_;
}


inline Foam::label Foam::ensightCells::total(const enum elemType what) const
{
    return sizes_[what];
}


inline Foam::label Foam::ensightCells::size(const enum elemType what) const
{
    return slices_[what].size();
}


inline Foam::label Foam::ensightCells::offset(const enum elemType what) const
{
    return slices_[what].start();
}


inline const Foam::labelUList Foam::ensightCells::cellIds
(
    const enum elemType what
) const
{
    return address_[slices_[what]];
}


inline const Foam::labelUList& Foam::ensightCells::cellIds() const
{
    return address_;
}


inline Foam::label Foam::ensightCells::operator[](const label i) const
{
    return address_[i];
}


// ************************************************************************* //
