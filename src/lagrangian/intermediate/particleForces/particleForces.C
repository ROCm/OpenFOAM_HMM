/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2011 OpenCFD Ltd.
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
#include "SRFModel.H"

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

Foam::particleForces::particleForces(const fvMesh& mesh)
:
    mesh_(mesh),
    dict_(dictionary::null),
    g_(vector::zero),
    gradUPtr_(NULL),
    HdotGradHInterPtr_(NULL),
    gravity_(false),
    virtualMass_(false),
    Cvm_(0.0),
    pressureGradient_(false),
    paramagnetic_(false),
    magneticSusceptibility_(0.0),
    refFrame_(rfInertial),
    UName_("undefined_UName"),
    HdotGradHName_("undefined_HdotGradHName")
{}


Foam::particleForces::particleForces
(
    const fvMesh& mesh,
    const dictionary& dict,
    const vector& g,
    const bool readFields
)
:
    mesh_(mesh),
    dict_(dict.subOrEmptyDict("particleForces")),
    g_(g),
    gradUPtr_(NULL),
    HdotGradHInterPtr_(NULL),
    gravity_(false),
    virtualMass_(false),
    Cvm_(0.0),
    pressureGradient_(false),
    paramagnetic_(false),
    magneticSusceptibility_(0.0),
    refFrame_(rfInertial),
    UName_(dict_.lookupOrDefault<word>("UName", "U")),
    HdotGradHName_(dict_.lookupOrDefault<word>("HdotGradHName", "HdotGradH"))
{
    if (readFields)
    {
        dict_.lookup("gravity") >> gravity_;
        dict_.lookup("virtualMass") >> virtualMass_;
        dict_.lookup("pressureGradient") >> pressureGradient_;
        dict_.lookup("paramagnetic") >> paramagnetic_;

        if (virtualMass_)
        {
            dict_.lookup("Cvm") >> Cvm_;
        }

        if (paramagnetic_)
        {
            dict_.lookup("magneticSusceptibility") >> magneticSusceptibility_;
        }
    }

    if (dict_.found("referenceFrame"))
    {
        word rf(dict_.lookup("referenceFrame"));

        if (rf == "SRF")
        {
            refFrame_ = rfSRF;
        }
        else if (rf != "inertial")
        {
            FatalErrorIn
            (
                "Foam::particleForces::particleForces"
                "("
                    "const fvMesh&, "
                    "const dictionary&, "
                    "const vector&, "
                    "const bool"
                ")"
            )
                << "Unknown referenceFrame, options are inertial and SRF."
                << abort(FatalError);
        }
    }
}


Foam::particleForces::particleForces(const particleForces& f)
:
    mesh_(f.mesh_),
    dict_(f.dict_),
    g_(f.g_),
    gradUPtr_(f.gradUPtr_),
    HdotGradHInterPtr_(f.HdotGradHInterPtr_),
    gravity_(f.gravity_),
    virtualMass_(f.virtualMass_),
    Cvm_(f.Cvm_),
    pressureGradient_(f.pressureGradient_),
    paramagnetic_(f.paramagnetic_),
    magneticSusceptibility_(f.magneticSusceptibility_),
    refFrame_(f.refFrame_),
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
            const volVectorField& HdotGradH =
                mesh_.lookupObject<volVectorField>(HdotGradHName_);

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
    const tetIndices& tetIs,
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
        accelTot += rhoc/rho*(U & gradU[tetIs.cell()]);
    }

    return accelTot;
}


Foam::vector Foam::particleForces::calcNonCoupled
(
    const vector& position,
    const tetIndices& tetIs,
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
           *HdotGradHInter.interpolate(position, tetIs);

        // force is:

        // 4.0
        // *constant::mathematical::pi
        // *constant::electromagnetic::mu0.value()
        // *pow3(d/2)
        // *magneticSusceptibility_/(magneticSusceptibility_ + 3)
        // *HdotGradH[cellI];

        // which is divided by mass (pi*d^3*rho/6) to produce
        // acceleration
    }

    if (refFrame_ == rfSRF)
    {
        const SRF::SRFModel& srf =
            mesh_.lookupObject<SRF::SRFModel>("SRFProperties");

        const vector& omega = srf.omega().value();
        const vector& axis = srf.axis();

        vector r = position - axis*(axis & position);

        // Coriolis and centrifugal acceleration terms
        accelTot += 2*(U ^ omega) + (omega ^ (r ^ omega));
    }

    return accelTot;
}


// ************************************************************************* //
