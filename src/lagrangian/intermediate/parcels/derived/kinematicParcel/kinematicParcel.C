/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

#include "kinematicParcel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(kinematicParcel, 0);
    defineParticleTypeNameAndDebug(kinematicParcel, 0);
    defineParcelTypeNameAndDebug(kinematicParcel, 0);
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kinematicParcel::kinematicParcel
(
    KinematicCloud<kinematicParcel>& owner,
    const label typeId,
    const vector& position,
    const label celli,
    const scalar d0,
    const vector& U0,
    const scalar nParticle0,
    const constantProperties& constProps
)
:
    KinematicParcel<kinematicParcel>
    (
        owner,
        typeId,
        position,
        celli,
        d0,
        U0,
        nParticle0,
        constProps
    )
{}


Foam::kinematicParcel::kinematicParcel
(
    const Cloud<kinematicParcel>& cloud,
    Istream& is,
    bool readFields
)
:
    KinematicParcel<kinematicParcel>(cloud, is, readFields)
{}


// * * * * * * * * * * * * * * * *  Destructors  * * * * * * * * * * * * * * //

Foam::kinematicParcel::~kinematicParcel()
{}


// ************************************************************************* //
