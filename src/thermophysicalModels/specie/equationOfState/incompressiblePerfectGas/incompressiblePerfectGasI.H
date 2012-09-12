/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

#include "incompressiblePerfectGas.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline Foam::incompressiblePerfectGas::incompressiblePerfectGas
(
    const specie& sp
)
:
    specie(sp)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::incompressiblePerfectGas::incompressiblePerfectGas
(
    const word& name,
    const incompressiblePerfectGas& pg
)
:
    specie(name, pg)
{}


inline Foam::autoPtr<Foam::incompressiblePerfectGas>
Foam::incompressiblePerfectGas::clone() const
{
    return autoPtr<incompressiblePerfectGas>
    (
        new incompressiblePerfectGas(*this)
    );
}


inline Foam::autoPtr<Foam::incompressiblePerfectGas>
Foam::incompressiblePerfectGas::New
(
    Istream& is
)
{
    return autoPtr<incompressiblePerfectGas>(new incompressiblePerfectGas(is));
}


inline Foam::autoPtr<Foam::incompressiblePerfectGas>
Foam::incompressiblePerfectGas::New
(
    const dictionary& dict
)
{
    return autoPtr<incompressiblePerfectGas>
    (
        new incompressiblePerfectGas(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::scalar Foam::incompressiblePerfectGas::rho
(
    scalar p,
    scalar T
) const
{
    return pRef_/(R()*T);
}


inline Foam::scalar Foam::incompressiblePerfectGas::psi(scalar, scalar T) const
{
    return 0.0;
}


inline Foam::scalar Foam::incompressiblePerfectGas::Z(scalar, scalar) const
{
    return 0.0;
}


inline Foam::scalar Foam::incompressiblePerfectGas::cpMcv(scalar, scalar) const
{
    return this->RR;
}



// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void Foam::incompressiblePerfectGas::operator+=
(
    const incompressiblePerfectGas& pg
)
{
    specie::operator+=(pg);
}


inline void Foam::incompressiblePerfectGas::operator-=
(
    const incompressiblePerfectGas& pg
)
{
    specie::operator-=(pg);
}


inline void Foam::incompressiblePerfectGas::operator*=(const scalar s)
{
    specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline Foam::incompressiblePerfectGas Foam::operator+
(
    const incompressiblePerfectGas& pg1,
    const incompressiblePerfectGas& pg2
)
{
    return incompressiblePerfectGas
    (
        static_cast<const specie&>(pg1)
      + static_cast<const specie&>(pg2)
    );
}


inline Foam::incompressiblePerfectGas Foam::operator-
(
    const incompressiblePerfectGas& pg1,
    const incompressiblePerfectGas& pg2
)
{
    return incompressiblePerfectGas
    (
        static_cast<const specie&>(pg1)
      - static_cast<const specie&>(pg2)
    );
}


inline Foam::incompressiblePerfectGas Foam::operator*
(
    const scalar s,
    const incompressiblePerfectGas& pg
)
{
    return incompressiblePerfectGas(s*static_cast<const specie&>(pg));
}


inline Foam::incompressiblePerfectGas Foam::operator==
(
    const incompressiblePerfectGas& pg1,
    const incompressiblePerfectGas& pg2
)
{
    return pg2 - pg1;
}


// ************************************************************************* //
