/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2005 OpenCFD Ltd.
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

inline const Foam::scalar Foam::pairPotential::m() const
{
    return m_;
}


inline const Foam::scalar Foam::pairPotential::gamma() const
{
    return gamma_;
}


inline const Foam::scalar Foam::pairPotential::rm() const
{
    return rm_;
}


inline const Foam::scalar Foam::pairPotential::epsilon() const
{
    return epsilon_;
}


inline const Foam::scalar Foam::pairPotential::rCut() const
{
    return rCut_;
}


inline const Foam::scalar Foam::pairPotential::rCutSqr() const
{
    return rCutSqr_;
}


inline const Foam::scalar Foam::pairPotential::rMin() const
{
    return rMin_;
}


inline const Foam::scalar Foam::pairPotential::dr() const
{
    return dr_;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline bool Foam::operator==
(
    const pairPotential& a,
    const pairPotential& b
)
{
    return
        (a.m() == b.m())
     && (a.gamma() == b.gamma())
     && (a.rm() == b.rm())
     && (a.epsilon() == b.epsilon())
     && (a.rCut() == b.rCut());
}


inline bool Foam::operator!=
(
    const pairPotential& a,
    const pairPotential& b
)
{
    return !(a == b);
}


// ************************************************************************* //
