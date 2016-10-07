/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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

#include "thermoCoupleProbes.H"
#include "mathematicalConstants.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(thermoCoupleProbes, 0);
    addToRunTimeSelectionTable(functionObject, thermoCoupleProbes, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::thermoCoupleProbes::readDict(const dictionary& dict)
{
    rho_ = readScalar(dict.lookup("rho"));
    Cp_ = readScalar(dict.lookup("Cp"));
    d_ = readScalar(dict.lookup("d"));
    epsilon_ = readScalar(dict.lookup("epsilon"));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::thermoCoupleProbes::thermoCoupleProbes
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    const bool loadFromFiles,
    const bool doFindElements
)
:
    probes(name, runTime, dict, loadFromFiles, doFindElements),
    ODESystem(),
    radName_(dict.lookup("radName")),
    thermo_
    (
        mesh_.lookupObject<fluidThermo>(basicThermo::dictName)
    ),
    odeSolver_(ODESolver::New(*this, dict)),
    Ttc_(this->size(), 0.0)
{

    readDict(dict);

    // Check if the property exist (resume old calculation)
    // or of it is new.
    if (foundProperty(typeName))
    {
        const dictionary& dict =
            this->stateDict().subDict(this->name()).subDict(typeName);

        dict.lookup("Tc") >> Ttc_;
    }
    else
    {
       Ttc_ = probes::sample(thermo_.T());
    }

    // Initialize thermocouple at initial T (room temperature)

    if (doFindElements)
    {
        // Find the elements
        findElements(mesh_);

        // Open the probe streams
        prepare();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::thermoCoupleProbes::~thermoCoupleProbes()
{}


Foam::label Foam::thermoCoupleProbes::nEqns() const
{
    return this->size();
}


void Foam::thermoCoupleProbes::derivatives
(
    const scalar x,
    const scalarField& y,
    scalarField& dydx
) const
{
    scalarField G(y.size(), 0.0);
    scalarField Tc(y.size(), 0.0);
    scalarField Uc(y.size(), 0.0);
    scalarField rhoc(y.size(), 0.0);
    scalarField muc(y.size(), 0.0);
    scalarField Cpc(y.size(), 0.0);
    scalarField kappac(y.size(), 0.0);

    if (radName_ != "none")
    {
        G = sample(mesh_.lookupObject<volScalarField>(radName_));
    }

    Tc = probes::sample(thermo_.T());

    Uc = mag(this->sample(mesh_.lookupObject<volVectorField>("U")));

    rhoc = this->sample(thermo_.rho()());
    kappac = this->sample(thermo_.kappa()());
    muc = this->sample(thermo_.mu()());
    Cpc = this->sample(thermo_.Cp()());

    scalarField Re(rhoc*Uc*d_/(muc + ROOTVSMALL));
    scalarField Pr(Cpc*muc/(kappac + ROOTVSMALL));
    //scalarField Nu(2.0 + 0.6*sqrt(Re)*cbrt(Pr));
    scalarField Nu(2.0 + (0.4*sqrt(Re) + 0.06*pow(Re, 2/3))*pow(Pr, 0.4));
    scalarField htc(Nu*kappac/d_);

    const scalar sigma = physicoChemical::sigma.value();

    scalar area = 4*constant::mathematical::pi*sqr(d_/2);
    scalar volume = (4/3)*constant::mathematical::pi*pow(d_/2, 3);

    dydx =
        (epsilon_*(G/4 - sigma*pow(y, 4.0))*area + htc*(Tc - y)*area)
      / (rho_*Cp_*volume);
}


void Foam::thermoCoupleProbes::jacobian
(
    const scalar x,
    const scalarField& y,
    scalarField& dfdt,
    scalarSquareMatrix& dfdy
) const
{
    derivatives(x, y, dfdt);
    for (label i=0; i<nEqns(); i++)
    {
        for (label j=0; j<nEqns(); j++)
        {
            dfdy(i, j) = 0.0;
        }
    }

}


bool Foam::thermoCoupleProbes::write()
{
    if (this->size())
    {
        sampleAndWrite<scalar>(thermo_.T());

        dictionary probeDict;
        probeDict.add("Tc", Ttc_);
        setProperty(typeName, probeDict);
        return true;
    }
     return false;
}


bool Foam::thermoCoupleProbes::execute()
{
    if (this->size())
    {
        scalar dt = mesh_.time().deltaTValue();
        scalar t = mesh_.time().value();
        odeSolver_->solve(t, t+dt, Ttc_, dt);
        return true;
    }

    return false;
}


bool Foam::thermoCoupleProbes::read(const dictionary& dict)
{
    readDict(dict);
    probes::read(dict);
    return true;
}


// ************************************************************************* //
