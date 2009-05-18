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

#include "COxidationMurphyShaddix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::label Foam::COxidationMurphyShaddix::maxIters_ = 1000;

Foam::scalar Foam::COxidationMurphyShaddix::tolerance_ = 1e-06;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::COxidationMurphyShaddix::COxidationMurphyShaddix
(
    const dictionary& dict,
    ReactingMultiphaseCloud<coalParcel>& owner
)
:
    SurfaceReactionModel<ReactingMultiphaseCloud<coalParcel> >
    (
        dict,
        owner,
        typeName
    ),
    D0_(dimensionedScalar(this->coeffDict().lookup("D0")).value()),
    rho0_(dimensionedScalar(this->coeffDict().lookup("rho0")).value()),
    T0_(dimensionedScalar(this->coeffDict().lookup("T0")).value()),
    Dn_(dimensionedScalar(this->coeffDict().lookup("Dn")).value()),
    A_(dimensionedScalar(this->coeffDict().lookup("A")).value()),
    E_(dimensionedScalar(this->coeffDict().lookup("E")).value()),
    n_(dimensionedScalar(this->coeffDict().lookup("n")).value()),
    WVol_(dimensionedScalar(this->coeffDict().lookup("WVol")).value()),
    dH0_(dimensionedScalar(this->coeffDict().lookup("dH0")).value()),
    TMaxdH_(dimensionedScalar(this->coeffDict().lookup("TmaxdH")).value()),
    CsLocalId_(-1),
    O2GlobalId_(-1),
    CO2GlobalId_(-1)
{
    // Determine carrier phase Ids of reactants/products
    label idGas = owner.composition().idGas();
    O2GlobalId_ = owner.composition().globalId(idGas, "O2");
    CO2GlobalId_ = owner.composition().globalId(idGas, "CO2");

    // Determine Cs ids
    label idSolid = owner.composition().idSolid();
    CsLocalId_ = owner.composition().localId(idSolid, "Cs");
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::COxidationMurphyShaddix::~COxidationMurphyShaddix()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::COxidationMurphyShaddix::active() const
{
    return true;
}


Foam::scalar Foam::COxidationMurphyShaddix::calculate
(
    const scalar dt,
    const label cellI,
    const scalar d,
    const scalar T,
    const scalar Tc,
    const scalar pc,
    const scalar rhoc,
    const scalar mass,
    const scalarField& YGas,
    const scalarField& YLiquid,
    const scalarField& YSolid,
    const scalarField& YMixture,
    const scalarField& dMassVolatile,
    scalarField& dMassGas,
    scalarField& dMassLiquid,
    scalarField& dMassSolid,
    scalarField& dMassSRCarrier
) const
{
    // Fraction of remaining combustible material
    const scalar fComb = YMixture[coalParcel::SLD]*YSolid[CsLocalId_];

    // Surface combustion until combustible fraction is consumed
    if (fComb < SMALL)
    {
        return 0.0;
    }

    // Cell carrier phase O2 species density [kg/m^3]
    const scalar rhoO2 =
        owner().carrierThermo().composition().Y(O2GlobalId_)[cellI]
       *rhoc;

    if (rhoO2 < SMALL)
    {
        return 0.0;
    }

    // Particle surface area [m^2]
    const scalar Ap = mathematicalConstant::pi*sqr(d);

//    // Film temperature [K]
//    const scalar Tf = (T + Tc)/2;

    // Calculate diffision constant at continuous phase temperature
    // and density [m^2/s]
    const scalar D = D0_*(rho0_/rhoc)*pow(Tc/T0_, Dn_);

    // Molecular weight of O2 [kg/kmol]
    const scalar WO2 = owner().composition().carrierSpecies()[O2GlobalId_].W();

    // Molecular weight of CO2 [kg/kmol]
    const scalar WCO2 =
        owner().composition().carrierSpecies()[CO2GlobalId_].W();

    // Molecular weight of C [kg/kmol]
    const scalar WC = WCO2 - WO2;

    // Far field partial pressure O2 [Pa]
    const scalar ppO2 = rhoO2/WO2*specie::RR*Tc;

    // Molar emission rate of volatiles per unit surface area
    const scalar qVol = sum(dMassVolatile)/(WVol_*Ap);

    // Total molar concentration of the carrier phase [kmol/m^3]
    const scalar C = pc/(specie::RR*Tc);

    if (debug)
    {
        Pout<< "mass  = " << mass << nl
            << "fComb = " << fComb << nl
            << "WC    = " << WC << nl
            << "Ap    = " << Ap << nl
            << "dt    = " << dt << nl
            << "C     = " << C << nl
            << endl;
    }

    // Molar reaction rate per unit surface area [kmol/(m^2.s)]
    scalar qCsOld = 0;
    scalar qCs = 1;

    const scalar qCsLim = mass*fComb/(WC*Ap*dt);

    if (debug)
    {
        Pout << "qCsLim = " << qCsLim << endl;
    }

    label iter = 0;
    while ((mag(qCs - qCsOld)/qCs > tolerance_) && (iter <= maxIters_))
    {
        qCsOld = qCs;
        const scalar PO2Surface = ppO2*exp(-(qCs + qVol)*d/(2*C*D));
        qCs = A_*exp(-E_/(specie::RR*T))*pow(PO2Surface, n_);
        qCs = (100.0*qCs + iter*qCsOld)/(100.0 + iter);
        qCs = min(qCs, qCsLim);

        if (debug)
        {
            Pout<< "iter = " << iter
                << ", qCsOld = " << qCsOld
                << ", qCs = " << qCs
                << nl << endl;
        }

        iter++;
    }

    if (iter > maxIters_)
    {
        WarningIn("scalar Foam::COxidationMurphyShaddix::calculate(...)")
            << "iter limit reached (" << maxIters_ << ")" << nl << endl;
    }

    // Calculate the number of molar units reacted
    scalar dOmega = qCs*Ap*dt;

    // Add to carrier phase mass transfer
    dMassSRCarrier[O2GlobalId_] += -dOmega*WO2;
    dMassSRCarrier[CO2GlobalId_] += dOmega*WCO2;

    // Add to particle mass transfer
    dMassSolid[CsLocalId_] = dOmega*WC;

    // Heat of reaction
    return dOmega*dH0_*cos(mathematicalConstant::pi/2*T/TMaxdH_);
}


// ************************************************************************* //
