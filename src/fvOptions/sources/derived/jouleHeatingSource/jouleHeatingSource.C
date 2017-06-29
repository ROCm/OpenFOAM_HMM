/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenCFD Ltd.
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

#include "jouleHeatingSource.H"
#include "fvMatrices.H"
#include "fvmLaplacian.H"
#include "fvcGrad.H"
#include "zeroGradientFvPatchField.H"
#include "basicThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(jouleHeatingSource, 0);

    addToRunTimeSelectionTable
    (
        option,
        jouleHeatingSource,
        dictionary
    );
}
}


const Foam::word Foam::fv::jouleHeatingSource::sigmaName(typeName + ":sigma");


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::coordinateSystem& Foam::fv::jouleHeatingSource::coordSys() const
{
    if (!coordSysPtr_.valid())
    {
        FatalErrorInFunction
            << "Co-ordinate system invalid"
            << abort(FatalError);
    }

    return coordSysPtr_();
}


Foam::tmp<Foam::volSymmTensorField>
Foam::fv::jouleHeatingSource::transformSigma
(
    const volVectorField& sigmaLocal
) const
{
    tmp<volSymmTensorField> tsigma
    (
        new volSymmTensorField
        (
            IOobject
            (
                sigmaName,
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedSymmTensor("0", sigmaLocal.dimensions(), Zero),
            zeroGradientFvPatchField<symmTensor>::typeName
        )
    );

    volSymmTensorField& sigma = tsigma.ref();
    sigma.primitiveFieldRef() = coordSys().R().transformVector(sigmaLocal);

    sigma.correctBoundaryConditions();

    return tsigma;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::jouleHeatingSource::jouleHeatingSource
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    option(sourceName, modelType, dict, mesh),
    TName_("T"),
    V_
    (
        IOobject
        (
            typeName + ":V",
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    ),
    anisotropicElectricalConductivity_(false),
    scalarSigmaVsTPtr_(nullptr),
    vectorSigmaVsTPtr_(nullptr),
    coordSysPtr_(nullptr),
    curTimeIndex_(-1)
{
    // Set the field name to that of the energy field from which the temperature
    // is obtained

    const basicThermo& thermo =
        mesh_.lookupObject<basicThermo>(basicThermo::dictName);

    fieldNames_.setSize(1, thermo.he().name());

    applied_.setSize(fieldNames_.size(), false);

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fv::jouleHeatingSource::~jouleHeatingSource()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

void Foam::fv::jouleHeatingSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    DebugInfo<< name() << ": applying source to " << eqn.psi().name() << endl;

    if (curTimeIndex_ != mesh_.time().timeIndex())
    {
        if (anisotropicElectricalConductivity_)
        {
            // Update sigma as a function of T if required
            const volVectorField& sigmaLocal = updateSigma(vectorSigmaVsTPtr_);

            tmp<volSymmTensorField> sigma = transformSigma(sigmaLocal);

            // Solve the electrical potential equation
            fvScalarMatrix VEqn(fvm::laplacian(sigma, V_));
            VEqn.relax();
            VEqn.solve();
        }
        else
        {
            // Update sigma as a function of T if required
            const volScalarField& sigma = updateSigma(scalarSigmaVsTPtr_);

            // Solve the electrical potential equation
            fvScalarMatrix VEqn(fvm::laplacian(sigma, V_));
            VEqn.relax();
            VEqn.solve();
        }

        curTimeIndex_ = mesh_.time().timeIndex();
    }

    // Add the Joule heating contribution

    const volVectorField gradV(fvc::grad(V_));
    if (anisotropicElectricalConductivity_)
    {
        const volVectorField& sigmaLocal =
            mesh_.lookupObject<volVectorField>(sigmaName);

        tmp<volSymmTensorField> sigma = transformSigma(sigmaLocal);

        eqn += (sigma & gradV) & gradV;
    }
    else
    {
        const volScalarField& sigma =
            mesh_.lookupObject<volScalarField>(sigmaName);

        eqn += (sigma*gradV) & gradV;
    }
}


// ************************************************************************* //
