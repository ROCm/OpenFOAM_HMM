/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "electricPotential.H"
#include "fvc.H"
#include "fvm.H"
#include "calculatedFvPatchField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(electricPotential, 0);
    addToRunTimeSelectionTable(functionObject, electricPotential, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::volScalarField&
Foam::functionObjects::electricPotential::operandField()
{
    if (!foundObject<volScalarField>(fieldName_))
    {
        auto tfldPtr = tmp<volScalarField>::New
        (
            IOobject
            (
                fieldName_,
                mesh_.time().timeName(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh_
        );
        store(fieldName_, tfldPtr);
    }

    return lookupObjectRef<volScalarField>(fieldName_);
}


Foam::tmp<Foam::volScalarField>
Foam::functionObjects::electricPotential::sigma() const
{
    const IOobject sigmaIO
    (
        IOobject::scopedName(typeName, "sigma"),
        mesh_.time().timeName(),
        mesh_.time(),
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
    );

    if (phases_.size())
    {
        tmp<volScalarField> tsigma = phases_[0]*sigmas_[0];

        for (label i = 1; i < phases_.size(); ++i)
        {
            tsigma.ref() += phases_[i]*sigmas_[i];
        }

        return tmp<volScalarField>::New
        (
            sigmaIO,
            tsigma,
            calculatedFvPatchField<scalar>::typeName
        );
    }

    return tmp<volScalarField>::New
    (
        sigmaIO,
        mesh_,
        sigma_,
        calculatedFvPatchField<scalar>::typeName
    );
}


Foam::tmp<Foam::volScalarField>
Foam::functionObjects::electricPotential::epsilonm() const
{
    // Vacuum permittivity (aka the electric constant) [A^2 s^4/(kg m^3)]
    const dimensionedScalar epsilon0
    (
        sqr(dimCurrent)*pow4(dimTime)/(dimMass*pow3(dimLength)),
        8.8541878128e-12    // CODATA value
    );

    const IOobject epsilonrIO
    (
        IOobject::scopedName(typeName, "epsilonr"),
        mesh_.time().timeName(),
        mesh_.time(),
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false
    );

    if (phases_.size())
    {
        tmp<volScalarField> tepsilonr = phases_[0]*epsilonrs_[0];

        for (label i = 1; i < phases_.size(); ++i)
        {
            tepsilonr.ref() += phases_[i]*epsilonrs_[i];
        }

        return tmp<volScalarField>::New
        (
            epsilonrIO,
            epsilon0*tepsilonr,
            calculatedFvPatchField<scalar>::typeName
        );
    }

    return tmp<volScalarField>::New
    (
        epsilonrIO,
        mesh_,
        epsilon0*epsilonr_,
        calculatedFvPatchField<scalar>::typeName
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::electricPotential::electricPotential
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    phasesDict_(dict.subOrEmptyDict("phases")),
    phaseNames_(),
    phases_(),
    sigmas_(),
    sigma_
    (
        dimensionedScalar
        (
            sqr(dimCurrent)*pow3(dimTime)/(dimMass*pow3(dimLength)),
            dict.getCheckOrDefault<scalar>
            (
                "sigma",
                scalar(1),
                scalarMinMax::ge(SMALL)
            )
        )
    ),
    epsilonrs_(),
    epsilonr_
    (
        dimensionedScalar
        (
            dimless,
            dict.getCheckOrDefault<scalar>
            (
                "epsilonr",
                scalar(1),
                scalarMinMax::ge(SMALL)
            )
        )
    ),
    fieldName_
    (
        dict.getOrDefault<word>
        (
            "field",
            IOobject::scopedName(typeName, "V")
        )
    ),
    nCorr_(1),
    writeDerivedFields_(false)
{
    read(dict);

    // Force creation of transported field so any BCs using it can
    // look it up
    volScalarField& eV = operandField();
    eV.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::electricPotential::read(const dictionary& dict)
{
    if (fvMeshFunctionObject::read(dict))
    {
        Log << type() << " read: " << name() << endl;

        dict.readIfPresent("sigma", sigma_);
        dict.readIfPresent("epsilonr", epsilonr_);
        dict.readIfPresent("nCorr", nCorr_);
        dict.readIfPresent("writeDerivedFields", writeDerivedFields_);

        // If flow is multiphase
        if (!phasesDict_.empty())
        {
            phaseNames_.setSize(phasesDict_.size());
            phases_.setSize(phasesDict_.size());
            sigmas_.setSize(phasesDict_.size());

            if (writeDerivedFields_)
            {
                epsilonrs_.setSize(phasesDict_.size());
            }

            label phasei = 0;
            for (const entry& dEntry : phasesDict_)
            {
                const word& key = dEntry.keyword();

                if (!dEntry.isDict())
                {
                    FatalIOErrorInFunction(phasesDict_)
                        << "Entry " << key << " is not a dictionary" << nl
                        << exit(FatalIOError);
                }

                const dictionary& subDict = dEntry.dict();

                phaseNames_[phasei] = key;

                sigmas_.set
                (
                    phasei,
                    new dimensionedScalar
                    (
                        sqr(dimCurrent)*pow3(dimTime)/(dimMass*pow3(dimLength)),
                        subDict.getCheck<scalar>
                        (
                            "sigma",
                            scalarMinMax::ge(SMALL)
                        )
                    )
                );

                if (writeDerivedFields_)
                {
                    epsilonrs_.set
                    (
                        phasei,
                        new dimensionedScalar
                        (
                            dimless,
                            subDict.getCheck<scalar>
                            (
                                "epsilonr",
                                scalarMinMax::ge(SMALL)
                            )
                        )
                    );
                }

                ++phasei;
            }

            forAll(phaseNames_, i)
            {
                phases_.set
                (
                    i,
                    mesh_.getObjectPtr<volScalarField>(phaseNames_[i])
                );
            }
        }

        return true;
    }

    return false;
}


bool Foam::functionObjects::electricPotential::execute()
{
    Log << type() << " execute: " << name() << endl;

    tmp<volScalarField> tsigma = this->sigma();
    const volScalarField& sigma = tsigma();

    volScalarField& eV = operandField();

    for (label i = 1; i <= nCorr_; ++i)
    {
        fvScalarMatrix eVEqn
        (
          - fvm::laplacian(sigma, eV)
        );

        eVEqn.relax();

        eVEqn.solve();
    }

    Log << endl;

    return true;
}


bool Foam::functionObjects::electricPotential::write()
{
    Log << type() << " write: " << name() << nl
        << tab << fieldName_
        << endl;

    volScalarField& eV = operandField();

    if (writeDerivedFields_)
    {
        // Write the electric field
        const volVectorField E
        (
            IOobject
            (
                IOobject::scopedName(typeName, "E"),
                mesh_.time().timeName(),
                mesh_.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            -fvc::grad(eV),
            calculatedFvPatchField<vector>::typeName
        );

        Log << tab << E.name() << endl;

        E.write();


        // Write the current density field
        tmp<volScalarField> tsigma = this->sigma();

        auto eJ = tmp<volVectorField>::New
        (
            IOobject
            (
                IOobject::scopedName(typeName, "J"),
                mesh_.time().timeName(),
                mesh_.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            -tsigma*fvc::grad(eV),
            calculatedFvPatchField<vector>::typeName
        );

        Log << tab << eJ().name() << endl;

        eJ->write();


        // Write the free-charge density field
        tmp<volScalarField> tepsilonm = this->epsilonm();

        auto erho = tmp<volScalarField>::New
        (
            IOobject
            (
                IOobject::scopedName(typeName, "rho"),
                mesh_.time().timeName(),
                mesh_.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            fvc::div(tepsilonm*E),
            calculatedFvPatchField<scalar>::typeName
        );

        Log << tab << erho().name() << endl;

        erho->write();
    }

    eV.write();

    return true;
}


// ************************************************************************* //
