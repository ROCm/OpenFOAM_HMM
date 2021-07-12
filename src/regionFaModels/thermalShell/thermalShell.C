/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "thermalShell.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFields.H"
#include "zeroGradientFaPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(thermalShell, 0);

addToRunTimeSelectionTable(thermalShellModel, thermalShell, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool thermalShell::init(const dictionary& dict)
{
    if (thickness_ > 0)
    {
        h_ = dimensionedScalar("thickness", dimLength, thickness_);
    }

    this->solution().readEntry("nNonOrthCorr", nNonOrthCorr_);

    return true;
}


tmp<areaScalarField> thermalShell::qr()
{
    IOobject io
    (
        "tqr",
        primaryMesh().time().timeName(),
        primaryMesh()
    );

    auto taqr =
        tmp<areaScalarField>::New
        (
            io,
            regionMesh(),
            dimensionedScalar(dimPower/dimArea, Zero)
        );

    if (qrName_ != "none")
    {
        auto& aqr = taqr.ref();

        const auto qr = primaryMesh().lookupObject<volScalarField>(qrName_);

        const volScalarField::Boundary& vqr = qr.boundaryField();

        aqr.primitiveFieldRef() = vsm().mapToSurface<scalar>(vqr);
    }

    return taqr;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void thermalShell::solveEnergy()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    const areaScalarField rhoCph(Cp()*rho()*h_);

    faScalarMatrix TEqn
    (
        fam::ddt(rhoCph, T_)
      - fam::laplacian(kappa()*h_, T_)
     ==
        qs_
      + qr()
      + faOptions()(h_, rhoCph, T_)
    );

    TEqn.relax();

    faOptions().constrain(TEqn);

    TEqn.solve();

    faOptions().correct(T_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

thermalShell::thermalShell
(
    const word& modelType,
    const fvPatch& patch,
    const dictionary& dict
)
:
    thermalShellModel(modelType, patch, dict),
    nNonOrthCorr_(1),
    thermo_(dict.subDict("thermo")),
    qs_
    (
        IOobject
        (
            "qs_" + regionName_,
            primaryMesh().time().timeName(),
            primaryMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimPower/dimArea, Zero)
    ),
    h_
    (
        IOobject
        (
            "h_" + regionName_,
            primaryMesh().time().timeName(),
            primaryMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh()
    ),
    qrName_(dict.getOrDefault<word>("qr", "none")),
    thickness_(dict.getOrDefault<scalar>("thickness", 0))
{
    init(dict);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void thermalShell::preEvolveRegion()
{}


void thermalShell::evolveRegion()
{
    nNonOrthCorr_ = solution().get<label>("nNonOrthCorr");

    for (int nonOrth = 0; nonOrth <= nNonOrthCorr_; ++nonOrth)
    {
        solveEnergy();
    }

    Info<< "T min/max   = " << min(T_) << ", " << max(T_) << endl;
}


const tmp<areaScalarField> thermalShell::Cp() const
{
    return tmp<areaScalarField>
    (
        new areaScalarField
        (
            IOobject
            (
                "Cps",
                primaryMesh().time().timeName(),
                primaryMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            regionMesh(),
            dimensionedScalar(dimEnergy/dimTemperature/dimMass, thermo_.Cp()),
            zeroGradientFaPatchScalarField::typeName
        )
    );
}


const tmp<areaScalarField> thermalShell::rho() const
{
    return tmp<areaScalarField>
    (
        new areaScalarField
        (
            IOobject
            (
                "rhos",
                primaryMesh().time().timeName(),
                primaryMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            regionMesh(),
            dimensionedScalar(dimDensity, thermo_.rho()),
            zeroGradientFaPatchScalarField::typeName
        )
    );
}


const tmp<areaScalarField> thermalShell::kappa() const
{
    return tmp<areaScalarField>
    (
        new areaScalarField
        (
            IOobject
            (
                "kappas",
                primaryMesh().time().timeName(),
                primaryMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            regionMesh(),
            dimensionedScalar
            (
                dimPower/dimLength/dimTemperature,
                thermo_.kappa()
            ),
            zeroGradientFaPatchScalarField::typeName
        )
    );
}


void thermalShell::info()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
