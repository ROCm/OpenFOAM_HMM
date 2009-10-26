/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

Foam::scalar Foam::uniformInterpolationTable::x0() const
{
    return x0_;
}


Foam::scalar Foam::uniformInterpolationTable::dx() const
{
    return dx_;
}


const Foam::Switch& Foam::uniformInterpolationTable::log10() const
{
    return log10_;
}


const Foam::Switch& Foam::uniformInterpolationTable::bound() const
{
    return bound_;
}


Foam::scalar& Foam::uniformInterpolationTable::x0()
{
    return x0_;
}


Foam::scalar& Foam::uniformInterpolationTable::dx()
{
    return dx_;
}


Foam::Switch& Foam::uniformInterpolationTable::log10()
{
    return log10_;
}


Foam::Switch& Foam::uniformInterpolationTable::bound()
{
    return bound_;
}


Foam::scalar Foam::uniformInterpolationTable::xMin() const
{
    return x0_;
}


Foam::scalar Foam::uniformInterpolationTable::xMax() const
{
    return x0_ + (size()-1)*dx_;
}


// ************************************************************************* //
