/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2010 OpenCFD Ltd.
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

#include "particleForces.H"
#include "fvMesh.H"
#include "volFields.H"
#include "fvcGrad.H"
#include "mathematicalConstants.H"
#include "electromagneticConstants.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::particleForces::deleteFields()
{
    if (gradUPtr_)
    {
        delete gradUPtr_;
        gradUPtr_ = NULL;
    }

    if (HdotGradHInterPtr_)
    {
        delete HdotGradHInterPtr_;
        HdotGradHInterPtr_ = NULL;
    }
}


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
    gradUPtr_(NULL),
    gravity_(dict_.lookup("gravity")),
    virtualMass_(dict_.lookup("virtualMass")),
    Cvm_(0.0),
    pressureGradient_(dict_.lookup("pressureGradient")),
    paramagnetic_(dict_.lookup("paramagnetic")),
    magneticSusceptibility_(0.0),
    UName_(dict_.lookupOrDefault<word>("U", "U")),
    HdotGradHName_(dict_.lookupOrDefault<word>("HdotGradH", "HdotGradH"))
{
    if (virtualMass_)
    {
        dict_.lookup("Cvm") >> Cvm_;
    }

    if (paramagnetic_)
    {
        dict_.lookup("magneticSusceptibility") >> magneticSusceptibility_;
    }
}


Foam::particleForces::particleForces(const particleForces& f)
:
    mesh_(f.mesh_),
    dict_(f.dict_),
    g_(f.g_),
    gradUPtr_(f.gradUPtr_),
    gravity_(f.gravity_),
    virtualMass_(f.virtualMass_),
    Cvm_(f.Cvm_),
    pressureGradient_(f.pressureGradient_),
    paramagnetic_(f.paramagnetic_),
    magneticSusceptibility_(f.magneticSusceptibility_),
    UName_(f.UName_),
    HdotGradHName_(f.HdotGradHName_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::particleForces::~particleForces()
{
    deleteFields();
}


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


Foam::scalar Foam::particleForces::Cvm() const
{
    return Cvm_;
}


Foam::Switch Foam::particleForces::pressureGradient() const
{
    return pressureGradient_;
}


Foam::Switch Foam::particleForces::paramagnetic() const
{
    return paramagnetic_;
}


Foam::scalar Foam::particleForces::magneticSusceptibility() const
{
    return magneticSusceptibility_;
}


const Foam::word& Foam::particleForces::UName() const
{
    return UName_;
}


const Foam::word& Foam::particleForces::HdotGradHName() const
{
    return HdotGradHName_;
}


void Foam::particleForces::cacheFields
(
    const bool store,
    const dictionary& interpolationSchemes
)
{
    if (store)
    {
        if (pressureGradient_)
        {
            const volVectorField& U =
                mesh_.lookupObject<volVectorField>(UName_);

            gradUPtr_ = fvc::grad(U).ptr();
        }

        if (paramagnetic_)
        {
            const volVectorField& HdotGradH = mesh_.lookupObject<volVectorField>
            (
                HdotGradHName_
            );

            HdotGradHInterPtr_ = interpolation<vector>::New
            (
                interpolationSchemes,
                HdotGradH
            ).ptr();
        }
    }
    else
    {
        deleteFields();
    }
}


Foam::vector Foam::particleForces::calcCoupled
(
    const vector& position,
    const label cellI,
    const scalar dt,
    const scalar rhoc,
    const scalar rho,
    const vector& Uc,
    const vector& U,
    const scalar d
) const
{
    vector accelTot = vector::zero;

    // Virtual mass force
    if (virtualMass_)
    {
        notImplemented
        (
            "Foam::particleForces::calcCoupled(...) - virtual mass force"
        );
//        accelTot += Cvm_*rhoc/rho*d(Uc - U)/dt;
    }

    // Pressure gradient force
    if (pressureGradient_)
    {
        const volTensorField& gradU = *gradUPtr_;
        accelTot += rhoc/rho*(U & gradU[cellI]);
    }

    return accelTot;
}


Foam::vector Foam::particleForces::calcNonCoupled
(
    const vector& position,
    const label cellI,
    const scalar dt,
    const scalar rhoc,
    const scalar rho,
    const vector& Uc,
    const vector& U,
    const scalar d
) const
{
    vector accelTot = vector::zero;

    // Gravity force
    if (gravity_)
    {
        accelTot += g_*(1.0 - rhoc/rho);
    }

    // Magnetic field force

    if (paramagnetic_)
    {
        const interpolation<vector>& HdotGradHInter = *HdotGradHInterPtr_;

        accelTot +=
            3.0*constant::electromagnetic::mu0.value()/rho
           *magneticSusceptibility_/(magneticSusceptibility_ + 3)
           *HdotGradHInter.interpolate(position, cellI);

        // force is:

        // 4.0
        // *constant::mathematical::pi
        // *constant::electromagnetic::mu0.value()
        // *pow3(d/2)
        // *magneticSusceptibility_/(magneticSusceptibility_ + 3)
        // *HdotGradH[cellI];

        // which is divided by mass ((4/3)*pi*r^3*rho) to produce
        // acceleration
    }

    return accelTot;
}


// ************************************************************************* //
