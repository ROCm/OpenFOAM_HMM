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

#include "ParticleTrackingData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class ParcelType>
Foam::label Foam::ParticleTrackingData<ParcelType>::PARTICLE_COUNT = 0;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::ParticleTrackingData<ParcelType>::ParticleTrackingData
(
    const Cloud<ParcelType>& cloud
)
:
    cloud_(cloud),
    origProc_(Pstream::myProcNo()),
    id_(PARTICLE_COUNT++)
{}


template<class ParcelType>
Foam::ParticleTrackingData<ParcelType>::ParticleTrackingData
(
    const ParticleTrackingData& ptd
)
:
    cloud_(ptd.cloud_),
    origProc_(ptd.origProc_),
    id_(ptd.id_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ParcelType>
Foam::ParticleTrackingData<ParcelType>::~ParticleTrackingData()
{}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

#include "ParticleTrackingDataIO.C"

// ************************************************************************* //
