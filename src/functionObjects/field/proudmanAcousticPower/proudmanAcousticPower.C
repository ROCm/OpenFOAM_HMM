/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "proudmanAcousticPower.H"
#include "volFields.H"
#include "basicThermo.H"
#include "turbulenceModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(proudmanAcousticPower, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        proudmanAcousticPower,
        dictionary
    );
}
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::functionObjects::proudmanAcousticPower::rhoScale
(
    const tmp<volScalarField>& fld
) const
{
    const basicThermo* thermoPtr =
        getObjectPtr<basicThermo>(basicThermo::dictName);

    if (thermoPtr)
    {
        return fld*thermoPtr->rho();
    }

    if (rhoInf_.value() < 0)
    {
        FatalErrorInFunction
            << type() << " " << name() << ": "
            << "Incompressible calculation assumed, but no reference density "
            << "set.  Please set the entry 'rhoInf' to an appropriate value"
            << nl
            << exit(FatalError);
    }

    return rhoInf_*fld;
}


Foam::tmp<Foam::volScalarField>
Foam::functionObjects::proudmanAcousticPower::a() const
{
    const basicThermo* thermoPtr =
        getObjectPtr<basicThermo>(basicThermo::dictName);

    if (thermoPtr)
    {
        const basicThermo& thermo = *thermoPtr;
        return sqrt(thermo.gamma()*thermo.p()/thermo.rho());
    }

    return
        tmp<volScalarField>::New
        (
            IOobject
            (
                scopedName("a"),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ
            ),
            mesh_,
            aRef_
        );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::proudmanAcousticPower::proudmanAcousticPower
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    alphaEps_(0.1),
    rhoInf_("0", dimDensity, -1),
    aRef_(dimVelocity, Zero)
{
    read(dict);

    volScalarField* PAPtr
    (
        new volScalarField
        (
            IOobject
            (
                scopedName("P_A"),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimPower/dimVolume, Zero)
        )
    );

    PAPtr->store();

    volScalarField* LPPtr
    (
        new volScalarField
        (
            IOobject
            (
                scopedName("L_P"),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimless, Zero)
        )
    );

    LPPtr->store();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::proudmanAcousticPower::read(const dictionary& dict)
{
    if (fvMeshFunctionObject::read(dict))
    {
        dict.readIfPresent("alphaEps", alphaEps_);
        rhoInf_.readIfPresent("rhoInf", dict);
        aRef_.readIfPresent("aRef", dict);

        return true;
    }

    return false;
}


bool Foam::functionObjects::proudmanAcousticPower::execute()
{
    const turbulenceModel& turb =
        lookupObject<turbulenceModel>(turbulenceModel::propertiesName);

    const volScalarField Mt(sqrt(2*turb.k())/a());

    volScalarField& P_A =
        mesh_.lookupObjectRef<volScalarField>(scopedName("P_A"));

    P_A = rhoScale(alphaEps_*turb.epsilon()*pow5(Mt));

    volScalarField& L_P =
        mesh_.lookupObjectRef<volScalarField>(scopedName("L_P"));

    L_P = 10*log10(P_A/dimensionedScalar("PRef", dimPower/dimVolume, 1e-12));

    return true;
}


bool Foam::functionObjects::proudmanAcousticPower::write()
{
    Log << type() << " " << name() << " write:" << nl;

    const volScalarField& P_A =
        mesh_.lookupObject<volScalarField>(scopedName("P_A"));

    Log << "    writing field " << P_A.name() << nl;

    P_A.write();

    const volScalarField& L_P =
        mesh_.lookupObject<volScalarField>(scopedName("L_P"));

    Log << "    writing field " << L_P.name() << nl;

    L_P.write();

    Log << endl;

    return true;
}


// ************************************************************************* //
