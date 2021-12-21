/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "KirchhoffShell.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFields.H"
#include "zeroGradientFaPatchFields.H"
#include "subCycle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(KirchhoffShell, 0);

addToRunTimeSelectionTable(vibrationShellModel, KirchhoffShell, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool KirchhoffShell::init(const dictionary& dict)
{
    this->solution().readEntry("nNonOrthCorr", nNonOrthCorr_);
    return true;
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void KirchhoffShell::solveDisplacement()
{
    if (debug)
    {
        InfoInFunction << endl;
    }

    const Time& time = primaryMesh().time();

    areaScalarField solidMass(rho()*h_);
    areaScalarField solidD(D()/solidMass);

    // Save old times
    areaScalarField w0(w_.oldTime());
    areaScalarField w00(w_.oldTime().oldTime());

    if (nSubCycles_ > 1)
    {
        // Restore the oldTime in sub-cycling
        w_.oldTime() = w0_;
        w_.oldTime().oldTime() = w00_;
        laplaceW_.oldTime() = laplaceW0_;
        laplace2W_.oldTime() = laplace2W0_;
     }

    for
    (
        subCycleTime wSubCycle
        (
            const_cast<Time&>(time),
            nSubCycles_
        );
       !(++wSubCycle).end();
    )
    {

        laplaceW_ = fac::laplacian(w_);
        laplace2W_ = fac::laplacian(laplaceW_);

        faScalarMatrix wEqn
        (
            fam::d2dt2(w_)
         +  f1_*fam::ddt(w_)
         -  f0_*sqrt(solidD)*fac::ddt(laplaceW_)
         +  solidD*(laplace2W_ + f2_*fac::ddt(laplace2W_))
        ==
            ps_/solidMass
          + faOptions()(solidMass, w_, dimLength/sqr(dimTime))
        );

        faOptions().constrain(wEqn);

        wEqn.solve();

        if (wSubCycle.index() >= wSubCycle.nSubCycles())
        {
            // Cache oldTimes inside the sub-cycling
            w0_ = w_.oldTime();
            w00_ = w_.oldTime().oldTime();
            laplaceW0_ = laplaceW_.oldTime();
            laplace2W0_ = laplace2W_.oldTime();

            // Update shell acceleration
            a_ = fac::d2dt2(w_);
        }
    }

    Info<< "ws_vibrationShell: "
        << "min = " << min(w_).value() << ", "
        << "max = " << max(w_).value() << endl;

    // Restore old time in main time
    w_.oldTime() = w0;
    w_.oldTime().oldTime() = w00;

    faOptions().correct(w_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

KirchhoffShell::KirchhoffShell
(
    const word& modelType,
    const fvPatch& patch,
    const dictionary& dict
)
:
    vibrationShellModel(modelType, patch, dict),
    f0_("f0", dimless, dict),
    f1_("f1", inv(dimTime), dict),
    f2_("f2", dimTime, dict),
    nNonOrthCorr_(1),
    nSubCycles_(1),
    ps_
    (
        IOobject
        (
            "ps_" + regionName_,
            primaryMesh().time().timeName(),
            primaryMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimPressure, Zero)
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
    laplaceW_
    (
        IOobject
        (
            "laplaceW_" + regionName_,
            primaryMesh().time().timeName(),
            primaryMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(inv(dimLength), Zero)
    ),
    laplace2W_
    (
        IOobject
        (
            "laplace2W_" + regionName_,
            primaryMesh().time().timeName(),
            primaryMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(inv(pow3(dimLength)), Zero)
    ),
    w0_
    (
        IOobject
        (
            "w0_" + regionName_,
            primaryMesh().time().timeName(),
            primaryMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimLength, Zero)
    ),
    w00_
    (
         IOobject
        (
            "w00_" + regionName_,
            primaryMesh().time().timeName(),
            primaryMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimLength, Zero)
    ),
    laplaceW0_
    (
         IOobject
        (
            "laplaceW0_" + regionName_,
            primaryMesh().time().timeName(),
            primaryMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(inv(dimLength), Zero)
    ),
    laplace2W0_
    (
         IOobject
        (
            "laplace2W0_" + regionName_,
            primaryMesh().time().timeName(),
            primaryMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(inv(pow3(dimLength)), Zero)
    )
{
    init(dict);
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void KirchhoffShell::preEvolveRegion()
{}


void KirchhoffShell::evolveRegion()
{
    nNonOrthCorr_ = solution().get<label>("nNonOrthCorr");
    nSubCycles_ = solution().get<label>("nSubCycles");

    for (int nonOrth=0; nonOrth<=nNonOrthCorr_; nonOrth++)
    {
        solveDisplacement();
    }
}


const tmp<areaScalarField> KirchhoffShell::D() const
{
    const dimensionedScalar E("E", dimForce/dimArea , solid().E());
    const dimensionedScalar nu("nu", dimless, solid().nu());

    return tmp<areaScalarField>(E*pow3(h_)/(12*(1 - sqr(nu))));
}


const tmp<areaScalarField> KirchhoffShell::rho() const
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
            dimensionedScalar("rho", dimDensity, solid().rho()),
            zeroGradientFaPatchScalarField::typeName
        )
    );
}

void KirchhoffShell::info()
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
