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

template<class Specie>
inline Foam::incompressiblePerfectGas<Specie>::incompressiblePerfectGas
(
    const Specie& sp
)
:
    Specie(sp)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::incompressiblePerfectGas<Specie>::incompressiblePerfectGas
(
    const word& name,
    const incompressiblePerfectGas<Specie>& pg
)
:
    Specie(name, pg)
{}


template<class Specie>
inline Foam::autoPtr<Foam::incompressiblePerfectGas<Specie> >
Foam::incompressiblePerfectGas<Specie>::clone() const
{
    return autoPtr<incompressiblePerfectGas<Specie> >
    (
        new incompressiblePerfectGas<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::incompressiblePerfectGas<Specie> >
Foam::incompressiblePerfectGas<Specie>::New
(
    Istream& is
)
{
    return autoPtr<incompressiblePerfectGas<Specie> >
    (
        new incompressiblePerfectGas<Specie>(is)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::incompressiblePerfectGas<Specie> >
Foam::incompressiblePerfectGas<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<incompressiblePerfectGas<Specie> >
    (
        new incompressiblePerfectGas<Specie>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::incompressiblePerfectGas<Specie>::rho
(
    scalar p,
    scalar T
) const
{
    return pRef_/(this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::incompressiblePerfectGas<Specie>::psi
(
    scalar,
    scalar T
) const
{
    return 0.0;
}


template<class Specie>
inline Foam::scalar Foam::incompressiblePerfectGas<Specie>::Z
(
    scalar,
    scalar
) const
{
    return 0.0;
}


template<class Specie>
inline Foam::scalar Foam::incompressiblePerfectGas<Specie>::cpMcv
(
    scalar,
    scalar
) const
{
    return this->RR;
}



// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::incompressiblePerfectGas<Specie>::operator+=
(
    const incompressiblePerfectGas<Specie>& pg
)
{
    Specie::operator+=(pg);
}


template<class Specie>
inline void Foam::incompressiblePerfectGas<Specie>::operator-=
(
    const incompressiblePerfectGas<Specie>& pg
)
{
    Specie::operator-=(pg);
}


template<class Specie>
inline void Foam::incompressiblePerfectGas<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::incompressiblePerfectGas<Specie> Foam::operator+
(
    const incompressiblePerfectGas<Specie>& pg1,
    const incompressiblePerfectGas<Specie>& pg2
)
{
    return incompressiblePerfectGas<Specie>
    (
        static_cast<const Specie&>(pg1)
      + static_cast<const Specie&>(pg2)
    );
}


template<class Specie>
inline Foam::incompressiblePerfectGas<Specie> Foam::operator-
(
    const incompressiblePerfectGas<Specie>& pg1,
    const incompressiblePerfectGas<Specie>& pg2
)
{
    return incompressiblePerfectGas<Specie>
    (
        static_cast<const Specie&>(pg1)
      - static_cast<const Specie&>(pg2)
    );
}


template<class Specie>
inline Foam::incompressiblePerfectGas<Specie> Foam::operator*
(
    const scalar s,
    const incompressiblePerfectGas<Specie>& pg
)
{
    return incompressiblePerfectGas<Specie>(s*static_cast<const Specie&>(pg));
}


template<class Specie>
inline Foam::incompressiblePerfectGas<Specie> Foam::operator==
(
    const incompressiblePerfectGas<Specie>& pg1,
    const incompressiblePerfectGas<Specie>& pg2
)
{
    return pg2 - pg1;
}


// ************************************************************************* //
