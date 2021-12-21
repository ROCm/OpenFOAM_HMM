/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 ENERCON GmbH
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

#include "atmPlantCanopyTurbSource.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(atmPlantCanopyTurbSource, 0);
    addToRunTimeSelectionTable(option, atmPlantCanopyTurbSource, dictionary);
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField::Internal>
Foam::fv::atmPlantCanopyTurbSource::calcPlantCanopyTerm
(
    const volVectorField::Internal& U
) const
{
    // (SP:Eq. 42)
    return 12.0*Foam::sqrt(Cmu_)*plantCd_()*leafAreaDensity_()*mag(U);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::atmPlantCanopyTurbSource::atmPlantCanopyTurbSource
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::cellSetOption(sourceName, modelType, dict, mesh),
    isEpsilon_(true),
    rhoName_(coeffs_.getOrDefault<word>("rho", "rho")),
    Cmu_(Zero),
    C1_(Zero),
    C2_(Zero),
    plantCd_
    (
        IOobject
        (
            "plantCd",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    leafAreaDensity_
    (
        IOobject
        (
            "leafAreaDensity",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    )
{
    const auto* turbPtr =
        mesh_.findObject<turbulenceModel>
        (
            turbulenceModel::propertiesName
        );

    if (!turbPtr)
    {
        FatalErrorInFunction
            << "Unable to find a turbulence model."
            << abort(FatalError);
    }

    fieldNames_.resize(1);

    tmp<volScalarField> tepsilon = turbPtr->epsilon();
    tmp<volScalarField> tomega = turbPtr->omega();

    if (!tepsilon.isTmp())
    {
        fieldNames_[0] = tepsilon().name();

        const dictionary& turbDict = turbPtr->coeffDict();
        Cmu_.read("Cmu", turbDict);
        C1_.read("C1", turbDict);
        C2_.read("C2", turbDict);
    }
    else if (!tomega.isTmp())
    {
        isEpsilon_ = false;
        fieldNames_[0] = tomega().name();

        const dictionary& turbDict = turbPtr->coeffDict();
        Cmu_.read("betaStar", turbDict);
    }
    else
    {
        FatalErrorInFunction
            << "Unable to find neither epsilon nor omega field." << nl
            << "atmPlantCanopyTurbSource needs either epsilon or omega field."
            << abort(FatalError);
    }

    fv::option::resetApplied();

    Log << "    Applying atmPlantCanopyTurbSource to: " << fieldNames_[0]
        << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::atmPlantCanopyTurbSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (isEpsilon_)
    {
        atmPlantCanopyTurbSourceEpsilon
        (
            geometricOneField(),
            geometricOneField(),
            eqn,
            fieldi
        );
    }
    else
    {
        atmPlantCanopyTurbSourceOmega
        (
            geometricOneField(),
            geometricOneField(),
            eqn,
            fieldi
        );
    }
}


void Foam::fv::atmPlantCanopyTurbSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (isEpsilon_)
    {
        atmPlantCanopyTurbSourceEpsilon(geometricOneField(), rho, eqn, fieldi);
    }
    else
    {
        atmPlantCanopyTurbSourceOmega(geometricOneField(), rho, eqn, fieldi);
    }
}


void Foam::fv::atmPlantCanopyTurbSource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (isEpsilon_)
    {
        atmPlantCanopyTurbSourceEpsilon(alpha, rho, eqn, fieldi);
    }
    else
    {
        atmPlantCanopyTurbSourceOmega(alpha, rho, eqn, fieldi);
    }
}


// ************************************************************************* //
