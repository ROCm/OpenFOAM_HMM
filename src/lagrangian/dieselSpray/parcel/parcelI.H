/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

inline const Foam::List<Foam::word>& Foam::parcel::liquidNames() const
{
    return liquidComponents_;
}


inline const Foam::List<Foam::word>& Foam::parcel::fuelNames() const
{
    return liquidComponents_;
}


inline Foam::scalar Foam::parcel::d() const
{
    return d_;
}


inline Foam::scalar& Foam::parcel::d()
{
    return d_;
}


inline Foam::scalar Foam::parcel::T() const
{
    return T_;
}


inline Foam::scalar& Foam::parcel::T()
{
    return T_;
}


inline Foam::scalar Foam::parcel::m() const
{
    return m_;
}


inline Foam::scalar& Foam::parcel::m()
{
    return m_;
}


inline Foam::scalar Foam::parcel::dev() const
{
    return y_;
}


inline Foam::scalar& Foam::parcel::dev()
{
    return y_;
}


inline Foam::scalar Foam::parcel::ddev() const
{
    return yDot_;
}


inline Foam::scalar& Foam::parcel::ddev()
{
    return yDot_;
}


inline Foam::scalar Foam::parcel::ct() const
{
    return ct_;
}


inline Foam::scalar& Foam::parcel::ct()
{
    return ct_;
}


inline Foam::scalar& Foam::parcel::ms()
{
    return ms_;
}


inline Foam::scalar Foam::parcel::ms() const
{
    return ms_;
}


inline Foam::scalar& Foam::parcel::tTurb()
{
    return tTurb_;
}


inline Foam::scalar Foam::parcel::tTurb() const
{
    return tTurb_;
}


inline Foam::scalar& Foam::parcel::liquidCore()
{
    return liquidCore_;
}


inline Foam::scalar Foam::parcel::liquidCore() const
{
    return liquidCore_;
}


inline Foam::scalar& Foam::parcel::injector()
{
    return injector_;
}


inline Foam::scalar Foam::parcel::injector() const
{
    return injector_;
}


inline const Foam::vector& Foam::parcel::U() const
{
    return U_;
}


inline Foam::vector& Foam::parcel::U()
{
    return U_;
}


inline const Foam::vector& Foam::parcel::Uturb() const
{
    return Uturb_;
}


inline Foam::vector& Foam::parcel::Uturb()
{
    return Uturb_;
}


inline const Foam::vector& Foam::parcel::n() const
{
    return n_;
}


inline Foam::vector& Foam::parcel::n()
{
    return n_;
}


inline const Foam::scalarField& Foam::parcel::X() const
{
    return X_;
}


inline Foam::scalarField& Foam::parcel::X()
{
    return X_;
}


inline Foam::scalar& Foam::parcel::tMom()
{
    return tMom_;
}


inline Foam::scalar Foam::parcel::tMom() const
{
    return tMom_;
}


inline Foam::vector Foam::parcel::Urel(const vector& U) const
{
    return U_ - U;
}


inline void Foam::parcel::correctNormal(const vector& sym)
{
    scalar d = position() & sym;
    n_ = position() - d*sym;
    n_ /= mag(n_);
}


// ************************************************************************* //
