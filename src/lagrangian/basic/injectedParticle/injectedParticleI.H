/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenCFD Ltd.
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::injectedParticle::injectedParticle
(
    const polyMesh& mesh,
    const vector& position,
    const label tag,
    const scalar soi,
    const scalar d,
    const vector& U,
    const bool doLocate
)
:
    particle(mesh, position, -1, -1, -1, doLocate),
    position_(position),
    tag_(tag),
    soi_(soi),
    d_(d),
    U_(U)
{}


inline Foam::label Foam::injectedParticle::tag() const
{
    return tag_;
}


inline Foam::scalar Foam::injectedParticle::soi() const
{
    return soi_;
}


inline Foam::scalar Foam::injectedParticle::d() const
{
    return d_;
}


inline const Foam::vector& Foam::injectedParticle::U() const
{
    return U_;
}


inline Foam::label& Foam::injectedParticle::tag()
{
    return tag_;
}


inline Foam::scalar& Foam::injectedParticle::soi()
{
    return soi_;
}


inline Foam::scalar& Foam::injectedParticle::d()
{
    return d_;
}


inline Foam::vector& Foam::injectedParticle::U()
{
    return U_;
}


// ************************************************************************* //
