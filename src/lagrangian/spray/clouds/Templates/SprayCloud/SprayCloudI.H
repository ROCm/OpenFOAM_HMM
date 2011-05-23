/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2011 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ParcelType>
inline const typename ParcelType::constantProperties&
Foam::SprayCloud<ParcelType>::constProps() const
{
    return constProps_;
}

template<class ParcelType>
inline const Foam::AtomizationModel<Foam::SprayCloud<ParcelType> >& Foam::SprayCloud<ParcelType>::atomization() const
{
    return atomizationModel_;
}

template<class ParcelType>
inline const Foam::BreakupModel<Foam::SprayCloud<ParcelType> >&
Foam::SprayCloud<ParcelType>::breakup() const
{
    return breakupModel_;
}

template<class ParcelType>
inline const Foam::CollisionModel<Foam::SprayCloud<ParcelType> >&
Foam::SprayCloud<ParcelType>::collision() const
{
    return collisionModel_;
}

template<class ParcelType>
inline const Foam::scalar& Foam::SprayCloud<ParcelType>::averageParcelMass() const
{
    return averageParcelMass_;
}

