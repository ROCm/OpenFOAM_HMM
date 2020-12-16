/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "LimitRange.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
inline Type Foam::Function1Types::LimitRange<Type>::value(const scalar t) const
{
    scalar tlim = min(max(t, min_), max_);

    return value_->value(tlim);
}


template<class Type>
Type Foam::Function1Types::LimitRange<Type>::integrate
(
    const scalar x1,
    const scalar x2
) const
{
    scalar xlim0 = min(max(x1, min_), max_);
    scalar xlim1 = min(max(x2, min_), max_);

    Type intValue = value_->integrate(xlim0, xlim1);

    if (x1 < min_)
    {
        intValue += (min(min_, x2) - x1)*this->value(min_);
    }

    if (x2 > max_)
    {
        intValue += (x2 - max(max_, x1))*this->value(max_);
    }

    return intValue;
}


// ************************************************************************* //
