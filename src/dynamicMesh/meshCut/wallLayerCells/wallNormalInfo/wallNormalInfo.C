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

#include "wallNormalInfo.H"

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Foam::Ostream& os,
    const Foam::wallNormalInfo& wDist
)
{
    if (os.format() == IOstream::ASCII)
    {
        os << wDist.normal();
    }
    else
    {
        os.write
        (
            reinterpret_cast<const char*>(&wDist.normal_),
            sizeof(vector)
        );
    }

    os.check(FUNCTION_NAME);
    return os;
}


Foam::Istream& Foam::operator>>(Foam::Istream& is, Foam::wallNormalInfo& wDist)
{
    if (is.format() == IOstream::ASCII)
    {
        is >> wDist.normal_;
    }
    else if (!is.checkLabelSize<>() || !is.checkScalarSize<>())
    {
        // Non-native label or scalar size
        is.beginRawRead();

        readRawScalar(is, wDist.normal_.data(), vector::nComponents);

        is.endRawRead();
    }
    else
    {
        is.read
        (
            reinterpret_cast<char*>(&wDist.normal_),
            sizeof(vector)
        );
    }

    is.check(FUNCTION_NAME);
    return is;
}


// ************************************************************************* //
