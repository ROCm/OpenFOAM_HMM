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

#include "particleForces.H"
#include "fvMesh.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::particleForces::particleForces
(
    const fvMesh& mesh,
    const dictionary& dict,
    const vector& g
)
:
    mesh_(mesh),
    dict_(dict.subDict("particleForces")),
    g_(g),
    gravity_(dict_.lookup("gravity")),
    virtualMass_(dict_.lookup("virtualMass")),
    Cvm_(0.0),
    pressureGradient_(dict_.lookup("pressureGradient")),
    gradUName_("unknown_gradUName")
{
    if (gravity_)
    {
        dict_.lookup("Cvm") >> Cvm_;
    }

    if (pressureGradient_)
    {
        dict_.lookup("gradU") >> gradUName_;
    }
}


Foam::particleForces::particleForces(const particleForces& f)
:
    mesh_(f.mesh_),
    dict_(f.dict_),
    g_(f.g_),
    gravity_(f.gravity_),
    virtualMass_(f.virtualMass_),
    Cvm_(f.Cvm_),
    pressureGradient_(f.pressureGradient_),
    gradUName_(f.gradUName_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::particleForces::~particleForces()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::dictionary& Foam::particleForces::dict() const
{
    return dict_;
}


const Foam::vector& Foam::particleForces::g() const
{
    return g_;
}


Foam::Switch Foam::particleForces::gravity() const
{
    return gravity_;
}


Foam::Switch Foam::particleForces::virtualMass() const
{
    return virtualMass_;
}


Foam::Switch Foam::particleForces::pressureGradient() const
{
    return pressureGradient_;
}


const Foam::word& Foam::particleForces::gradUName() const
{
    return gradUName_;
}


Foam::vector Foam::particleForces::calc
(
    const label cellI,
    const scalar dt,
    const scalar rhoc,
    const scalar rho,
    const vector& Uc,
    const vector& U
) const
{
    vector Ftot = vector::zero;

    // Gravity force
    if (gravity_)
    {
        Ftot += g_*(1.0 - rhoc/rho);
    }

    // Virtual mass force
    if (virtualMass_)
    {
        notImplemented("Foam::particleForces::calc(...) - virtualMass force");
//        Ftot += Cvm_*rhoc/rho*d(Uc - U)/dt;
    }

    // Pressure gradient force
    if (pressureGradient_)
    {
        const volSymmTensorField& gradU =
            mesh_.lookupObject<volSymmTensorField>(gradUName_);
        Ftot += rhoc/rho*(U & gradU[cellI]);
    }

    return Ftot;
}


// ************************************************************************* //

