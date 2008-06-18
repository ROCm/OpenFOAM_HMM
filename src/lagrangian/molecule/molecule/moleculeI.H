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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::molecule::molecule
(
    const Cloud<molecule>& c,
    const vector& position,
    const label celli,
    const scalar mass,
    const vector& U,
    const vector& A,
    const vector& tetherPosition,
    const label tethered,
    const label id
)
:
    Particle<molecule>(c, position, celli),
    mass_(mass),
    U_(U),
    A_(A),
    tetherPosition_(tetherPosition),
    potentialEnergy_(0.0),
    rf_(tensor::zero),
    tethered_(tethered),
    id_(id)
{}


inline Foam::molecule::trackData::trackData
(
    moleculeCloud& molCloud,
    label part
)
:
    Particle<molecule>::trackData(refCast<Cloud<molecule> >(molCloud)),
    molCloud_(molCloud),
    part_(part)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::molecule::id() const
{
    return id_;
}


inline Foam::scalar Foam::molecule::mass() const
{
    return mass_;
}


inline const Foam::vector& Foam::molecule::U() const
{
    return U_;
}


inline Foam::vector& Foam::molecule::U()
{
    return U_;
}


inline const Foam::vector& Foam::molecule::A() const
{
    return A_;
}


inline Foam::vector& Foam::molecule::A()
{
    return A_;
}


inline Foam::scalar Foam::molecule::potentialEnergy() const
{
    return potentialEnergy_;
}


inline Foam::scalar& Foam::molecule::potentialEnergy()
{
    return potentialEnergy_;
}


inline const Foam::tensor& Foam::molecule::rf() const
{
    return rf_;
}


inline Foam::tensor& Foam::molecule::rf()
{
    return rf_;
}


inline Foam::label Foam::molecule::tethered() const
{
    return tethered_;
}


inline const Foam::vector& Foam::molecule::tetherPosition() const
{
    return tetherPosition_;
}


inline Foam::moleculeCloud& Foam::molecule::trackData::molCloud()
{
    return molCloud_;
}


inline Foam::label Foam::molecule::trackData::part() const
{
    return part_;
}


// ************************************************************************* //
