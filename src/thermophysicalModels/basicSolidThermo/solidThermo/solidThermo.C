/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "solidThermo.H"
#include "fvMesh.H"
#include "volFields.H"
#include "HashTable.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(solidThermo, 0);
    defineRunTimeSelectionTable(solidThermo, mesh);
    defineRunTimeSelectionTable(solidThermo, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidThermo::solidThermo(const fvMesh& mesh)
:
    basicThermo(mesh),
    rho_
    (
        IOobject
        (
            "rho",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimMass/dimVolume
    )
{}


Foam::solidThermo::solidThermo
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    basicThermo(mesh, dict),
    rho_
    (
        IOobject
        (
            "rho",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimMass/dimVolume
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidThermo::~solidThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::volScalarField& Foam::solidThermo::T()
{
    return this->T_;
}


const Foam::volScalarField& Foam::solidThermo::T() const
{
    return this->T_;
}


const Foam::volScalarField& Foam::solidThermo::rhos() const
{
    return rho_;
}


Foam::volScalarField& Foam::solidThermo::rhos()
{
    return rho_;
}


const Foam::volScalarField& Foam::solidThermo::p() const
{
    return this->p_;
}


Foam::volScalarField& Foam::solidThermo::p()
{
    return this->p_;
}


const Foam::volScalarField& Foam::solidThermo::alpha() const
{
    return this->alpha_;
}


Foam::tmp<Foam::volScalarField> Foam::solidThermo::rho() const
{
    return tmp<volScalarField>(rho_);
}


bool Foam::solidThermo::read()
{
    return regIOobject::read();
}

// ************************************************************************* //
