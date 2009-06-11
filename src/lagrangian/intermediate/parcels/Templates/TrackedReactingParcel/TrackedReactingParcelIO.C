/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2009 OpenCFD Ltd.
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

#include "TrackedReactingParcel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template <class ParcelType>
Foam::TrackedReactingParcel<ParcelType>::TrackedReactingParcel
(
    const Cloud<ParcelType>& cloud,
    Istream& is,
    bool readFields
)
:
    ReactingParcel<ParcelType>(cloud, is, readFields),
    ParticleTrackingData<ParcelType>(cloud, is, readFields)
{}


template<class ParcelType>
void Foam::TrackedReactingParcel<ParcelType>::readFields
(
    ReactingCloud<ParcelType>& c
)
{
    if (!c.size())
    {
        return;
    }

    ReactingParcel<ParcelType>::readFields(c);
    ParticleTrackingData<ParcelType>::readFields(c);
}


template<class ParcelType>
void Foam::TrackedReactingParcel<ParcelType>::writeFields
(
    const ReactingCloud<ParcelType>& c
)
{
    ReactingParcel<ParcelType>::writeFields(c);
    ParticleTrackingData<ParcelType>::writeFields(c);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ParcelType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const TrackedReactingParcel<ParcelType>& p
)
{
    if (os.format() == IOstream::ASCII)
    {
        os  << static_cast<const ReactingParcel<ParcelType>&>(p)
            << static_cast<const ParticleTrackingData<ParcelType>&>(p);
    }
    else
    {
        os  << static_cast<const ReactingParcel<ParcelType>&>(p)
            << static_cast<const ParticleTrackingData<ParcelType>&>(p);
    }

    // Check state of Ostream
    os.check
    (
        "Ostream& operator<<"
        "("
            "Ostream&, "
            "const TrackedReactingParcel<ParcelType>&"
        ")"
    );

    return os;
}


// ************************************************************************* //
