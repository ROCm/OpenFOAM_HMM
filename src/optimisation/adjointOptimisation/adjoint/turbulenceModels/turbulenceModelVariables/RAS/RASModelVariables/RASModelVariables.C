/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2007-2019 PCOpt/NTUA
                            | Copyright (C) 2013-2019 FOSS GP
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

#include "RASModelVariables.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(RASModelVariables, 0);
defineRunTimeSelectionTable(RASModelVariables, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void RASModelVariables::allocateInitValues()
{
    if (solverControl_.storeInitValues())
    {
        Info<< "Storing initial values of turbulence variables" << endl;
        if (hasTMVar1_)
        {
            TMVar1InitPtr_.reset
            (
                new volScalarField
                (
                    TMVar1Inst().name()+"Init",TMVar1Inst()
                )
            );
        }

        if (hasTMVar2_)
        {
            TMVar2InitPtr_.reset
            (
                new volScalarField
                (
                    TMVar2Inst().name()+"Init",TMVar2Inst()
                )
            );
        }

        if (hasNut_)
        {
            nutInitPtr_.reset
            (
                new volScalarField
                (
                    nutRefInst().name()+"Init",nutRefInst()
                )
            );
        }
    }
}


void RASModelVariables::allocateMeanFields()
{
    if (solverControl_.average())
    {
        Info<< "Allocating mean values of turbulence variables" << endl;
        if (hasTMVar1_)
        {
            TMVar1MeanPtr_.reset
            (
                new volScalarField
                (
                    IOobject
                    (
                        TMVar1Inst().name()+"Mean",
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE
                    ),
                    TMVar1Inst()
                )
            );
        }
        if (hasTMVar2_)
        {
            TMVar2MeanPtr_.reset
            (
                new volScalarField
                (
                    IOobject
                    (
                        TMVar2Inst().name()+"Mean",
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE
                    ),
                    TMVar2Inst()
                )
            );
        }

        if (hasNut_)
        {
            nutMeanPtr_.reset
            (
                new volScalarField
                (
                    IOobject
                    (
                        nutRefInst().name()+"Mean",
                        mesh_.time().timeName(),
                        mesh_,
                        IOobject::READ_IF_PRESENT,
                        IOobject::AUTO_WRITE
                    ),
                    nutRefInst()
                )
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

RASModelVariables::RASModelVariables
(
    const fvMesh& mesh,
    const solverControl& SolverControl
)
:
    mesh_(mesh),
    solverControl_(SolverControl),
    hasTMVar1_(false),
    hasTMVar2_(false),
    hasNut_(false),
    hasDist_(false),
    TMVar1Ptr_(nullptr),
    TMVar2Ptr_(nullptr),
    nutPtr_(nullptr),
    dPtr_(nullptr),
    TMVar1BaseName_(word::null),
    TMVar2BaseName_(word::null),
    nutBaseName_("nut"),
    TMVar1InitPtr_(nullptr),
    TMVar2InitPtr_(nullptr),
    nutInitPtr_(nullptr),
    TMVar1MeanPtr_(nullptr),
    TMVar2MeanPtr_(nullptr),
    nutMeanPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<RASModelVariables> RASModelVariables::New
(
    const fvMesh& mesh,
    const solverControl& SolverControl
)
{
    // Get model name, but do not register the dictionary
    // otherwise it is registered in the database twice
    const word modelType
    (
        IOdictionary
        (
            IOobject
            (
                turbulenceModel::propertiesName,
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        ).subOrEmptyDict("RAS").lookupOrDefault<word>("RASModel", "laminar")
    );

    Info<< "Creating references for RASModel variables : " << modelType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown RASModelVariables type " << modelType << nl << nl
            << "Valid RASModelVariables types are :" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<RASModelVariables>(cstrIter()(mesh, SolverControl));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool RASModelVariables::hasTMVar1() const
{
    return hasTMVar1_;
}


bool RASModelVariables::hasTMVar2() const
{
    return hasTMVar2_;
}


bool RASModelVariables::hasNut() const
{
    return hasNut_;
}


bool RASModelVariables::hasDist() const
{
    return hasDist_;
}


const word& RASModelVariables::TMVar1BaseName() const
{
    return TMVar1BaseName_;
}


const word& RASModelVariables::TMVar2BaseName() const
{
    return TMVar2BaseName_;
}


const word& RASModelVariables::nutBaseName() const
{
    return nutBaseName_;
}


const volScalarField& RASModelVariables::TMVar1() const
{
    if (solverControl_.useAveragedFields())
    {
        return TMVar1MeanPtr_();
    }
    else
    {
        return *TMVar1Ptr_;
    }
}


volScalarField& RASModelVariables::TMVar1()
{
    if (solverControl_.useAveragedFields())
    {
        return TMVar1MeanPtr_();
    }
    else
    {
        return *TMVar1Ptr_;
    }
}


const volScalarField& RASModelVariables::TMVar2() const
{
    if (solverControl_.useAveragedFields())
    {
        return TMVar2MeanPtr_();
    }
    else
    {
        return *TMVar2Ptr_;
    }
}

volScalarField& RASModelVariables::TMVar2()
{
    if (solverControl_.useAveragedFields())
    {
        return TMVar2MeanPtr_();
    }
    else
    {
        return *TMVar2Ptr_;
    }
}

const volScalarField& RASModelVariables::nutRef() const
{
    if (solverControl_.useAveragedFields() && hasNut_)
    {
        return nutMeanPtr_();
    }
    else
    {
        return *nutPtr_;
    }
}


volScalarField& RASModelVariables::nutRef()
{
    if (solverControl_.useAveragedFields() && hasNut_)
    {
        return  nutMeanPtr_();
    }
    else
    {
        return *nutPtr_;
    }
}


const volScalarField& RASModelVariables::d() const
{
    return *dPtr_;
}


volScalarField& RASModelVariables::d()
{
    return *dPtr_;
}


const volScalarField& RASModelVariables::TMVar1Inst() const
{
    return *TMVar1Ptr_;
}


volScalarField& RASModelVariables::TMVar1Inst()
{
    return *TMVar1Ptr_;
}


const volScalarField& RASModelVariables::TMVar2Inst() const
{
    return *TMVar2Ptr_;
}


volScalarField& RASModelVariables::TMVar2Inst()
{
    return *TMVar2Ptr_;
}


const volScalarField& RASModelVariables::nutRefInst() const
{
    return *nutPtr_;
}


volScalarField& RASModelVariables::nutRefInst()
{
    return *nutPtr_;
}


tmp<volScalarField> RASModelVariables::nutJacobianVar1
(
    const singlePhaseTransportModel& laminarTransport
) const
{
    WarningInFunction
        << "jutJacobianVar1 not implemented for the current turbulence model."
        << "Returning zero field" << endl;

    tmp<volScalarField> nutJacobian
    (
        new volScalarField
        (
            IOobject
            (
                "nutJacobianVar1",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimless, Zero)
        )
    );
    return nutJacobian;
}


tmp<volScalarField> RASModelVariables::nutJacobianVar2
(
    const singlePhaseTransportModel& laminarTransport
) const
{
    WarningInFunction
        << "nutJacobianVar2 not implemented for the current turbulence model."
        << "Returning zero field" << endl;

    tmp<volScalarField> nutJacobian
    (
        new volScalarField
        (
            IOobject
            (
                "nutJacobianVar2",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimless, Zero)
        )
    );
    return nutJacobian;
}

void RASModelVariables::restoreInitValues()
{
    if (solverControl_.storeInitValues())
    {
        if (hasTMVar1_)
        {
            TMVar1Inst() == TMVar1InitPtr_();
        }
        if (hasTMVar2_)
        {
            TMVar2Inst() == TMVar2InitPtr_();
        }
        if (hasNut_)
        {
            nutRefInst() == nutInitPtr_();
        }
    }
}


void RASModelVariables::resetMeanFields()
{
    if (solverControl_.average())
    {
        Info<< "Reseting mean turbulent fields to zero" << endl;

        // Reset fields to zero
        if (hasTMVar1_)
        {
            TMVar1MeanPtr_() ==
                dimensionedScalar(TMVar1Inst().dimensions(), Zero);
        }
        if (hasTMVar2_)
        {
            TMVar2MeanPtr_() ==
                dimensionedScalar(TMVar2Inst().dimensions(), Zero);
        }
        if (hasNut_)
        {
            nutMeanPtr_() == dimensionedScalar(nutRefInst().dimensions(), Zero);
        }
    }
}


void RASModelVariables::computeMeanFields()
{
    if (solverControl_.doAverageIter())
    {
        const label iAverageIter = solverControl_.averageIter();
        scalar avIter(iAverageIter);
        scalar oneOverItP1 = 1./(avIter + 1);
        scalar mult = avIter*oneOverItP1;
        if (hasTMVar1_)
        {
            TMVar1MeanPtr_() ==
                TMVar1MeanPtr_()*mult + TMVar1Inst()*oneOverItP1;
        }
        if (hasTMVar2_)
        {
            TMVar2MeanPtr_() ==
                TMVar2MeanPtr_()*mult + TMVar2Inst()*oneOverItP1;
        }
        if (hasNut_)
        {
            nutMeanPtr_() == nutMeanPtr_()*mult + nutRefInst()*oneOverItP1;
        }
    }
}


tmp<volSymmTensorField> RASModelVariables::devReff
(
    const singlePhaseTransportModel& laminarTransport,
    const volVectorField& U
) const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -(laminarTransport.nu() + nutRef())*dev(twoSymm(fvc::grad(U)))
        )
    );
}


void RASModelVariables::correctBoundaryConditions
(
    const incompressible::turbulenceModel& turbulence
)
{
    if (hasTMVar1())
    {
        TMVar1Ptr_->correctBoundaryConditions();
        if (solverControl_.average())
        {
            TMVar1MeanPtr_().correctBoundaryConditions();
        }
    }

    if (hasTMVar2())
    {
        TMVar2Ptr_->correctBoundaryConditions();
        if (solverControl_.average())
        {
            TMVar2MeanPtr_().correctBoundaryConditions();
        }
    }

    if (hasNut())
    {
        nutPtr_->correctBoundaryConditions();
        if (solverControl_.average())
        {
            nutMeanPtr_().correctBoundaryConditions();
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
