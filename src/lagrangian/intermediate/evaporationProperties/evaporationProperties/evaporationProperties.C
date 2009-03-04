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

#include "evaporationProperties.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::evaporationProperties::evaporationProperties()
:
    name_("unknownSpecie"),
    Dab_(0.0),
    TvsPSat_()
{}


Foam::evaporationProperties::evaporationProperties
(
    const evaporationProperties& pp
)
:
    name_(pp.name_),
    Dab_(pp.Dab_),
    TvsPSat_(pp.TvsPSat_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::evaporationProperties::~evaporationProperties()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::word& Foam::evaporationProperties::name() const
{
    return name_;
}


Foam::scalar Foam::evaporationProperties::Dab() const
{
    return Dab_;
}


const Foam::DataEntry<Foam::scalar>&
Foam::evaporationProperties::TvsPSat() const
{
    return TvsPSat_();
}


// ************************************************************************* //

