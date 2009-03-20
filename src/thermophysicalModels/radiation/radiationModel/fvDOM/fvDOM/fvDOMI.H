/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

inline const Foam::radiation::radiativeIntensityRay&
Foam::radiation::fvDOM::RadIntRay(label i) const
{
    return  RadIntRay_[i];
}


inline const Foam::volScalarField& Foam::radiation::fvDOM::RadIntRayiLambdaj
(
    const label i,
    const label j
) const
{
    return  RadIntRay_[i].Ilambdaj(j);
}


inline Foam::label Foam::radiation::fvDOM::Ntheta() const
{
    return Ntheta_;
}


inline Foam::label Foam::radiation::fvDOM::Nphi() const
{
    return Nphi_;
}


inline Foam::label Foam::radiation::fvDOM::Ni() const
{
    return Ni_;
}


inline Foam::label Foam::radiation::fvDOM::lambdaj() const
{
    return lambdaj_;
}


inline const Foam::volScalarField& Foam::radiation::fvDOM::a() const
{
    return a_;
}


inline const Foam::volScalarField& Foam::radiation::fvDOM::aj
(
    const label i
) const
{
    return aj_[i];
}


inline const Foam::volScalarField& Foam::radiation::fvDOM::G() const
{
    return G_;
}


inline const Foam::volScalarField& Foam::radiation::fvDOM::Qr() const
{
    return Qr_;
}


virtual const Foam::radiation::blackBodyEmission&
Foam::radiation::fvDOM::blackBody() const
{
    return blackBody_;
}


// ************************************************************************* //
