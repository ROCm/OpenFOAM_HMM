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

#include "RASModel.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(RASModel, 0);
defineRunTimeSelectionTable(RASModel, dictionary);
addToRunTimeSelectionTable(turbulenceModel, RASModel, turbulenceModel);

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void RASModel::printCoeffs()
{
    if (printCoeffs_)
    {
        Info<< type() << "Coeffs" << coeffDict_ << nl
            << "wallFunctionCoeffs" << wallFunctionDict_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

RASModel::RASModel
(
    const word& type,
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& lamTransportModel
)
:
    turbulenceModel(U, phi, lamTransportModel),

    IOdictionary
    (
        IOobject
        (
            "RASProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    turbulence_(lookup("turbulence")),
    printCoeffs_(lookupOrDefault<Switch>("printCoeffs", false)),
    coeffDict_(subDict(type + "Coeffs")),

    wallFunctionDict_(subDict("wallFunctionCoeffs")),
    kappa_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kappa",
            wallFunctionDict_,
            0.4187
        )
    ),
    E_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "E",
            wallFunctionDict_,
            9.0
        )
    ),
    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            wallFunctionDict_,
            0.09
        )
    ),

    yPlusLam_(yPlusLam(kappa_.value(), E_.value())),

    k0_("k0", dimVelocity*dimVelocity, SMALL),
    epsilon0_("epsilon", k0_.dimensions()/dimTime, SMALL),
    epsilonSmall_("epsilonSmall", epsilon0_.dimensions(), SMALL),

    y_(mesh_)
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<RASModel> RASModel::New
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport
)
{
    word RASModelTypeName;

    // Enclose the creation of the turbulencePropertiesDict to ensure it is
    // deleted before the RASModel is created otherwise the dictionary
    // is entered in the database twice
    {
        IOdictionary turbulencePropertiesDict
        (
            IOobject
            (
                "RASProperties",
                U.time().constant(),
                U.db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        turbulencePropertiesDict.lookup("RASModel")
            >> RASModelTypeName;
    }

    Info<< "Selecting RAS turbulence model " << RASModelTypeName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(RASModelTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "RASModel::New(const volVectorField&, "
            "const surfaceScalarField&, transportModel&)"
        )   << "Unknown RASModel type " << RASModelTypeName
            << endl << endl
            << "Valid RASModel types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<RASModel>(cstrIter()(U, phi, transport));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

RASModel::~RASModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar RASModel::yPlusLam(const scalar kappa, const scalar E) const
{
    scalar ypl = 11.0;

    for (int i=0; i<10; i++)
    {
        ypl = log(max(E*ypl, 1))/kappa;
    }

    return ypl;
}


tmp<scalarField> RASModel::yPlus(const label patchNo) const
{
    const fvPatch& curPatch = mesh_.boundary()[patchNo];

    tmp<scalarField> tYp(new scalarField(curPatch.size()));
    scalarField& Yp = tYp();

    if (isType<wallFvPatch>(curPatch))
    {
        Yp = pow(Cmu_.value(), 0.25)
            *y_[patchNo]
            *sqrt(k()().boundaryField()[patchNo].patchInternalField())
            /nu().boundaryField()[patchNo];
    }
    else
    {
        WarningIn
        (
            "tmp<scalarField> RASModel::yPlus(const label patchNo) const"
        )   << "Patch " << patchNo << " is not a wall. Returning null field"
            << nl << endl;

        Yp.setSize(0);
    }

    return tYp;
}


void RASModel::correct()
{
    turbulenceModel::correct();

    if (turbulence_ && mesh_.changing())
    {
        y_.correct();
    }
}


bool RASModel::read()
{
    if (regIOobject::read())
    {
        lookup("turbulence") >> turbulence_;
        coeffDict_ = subDict(type() + "Coeffs");

        wallFunctionDict_ = subDict("wallFunctionCoeffs");
        kappa_.readIfPresent(wallFunctionDict_);
        E_.readIfPresent(wallFunctionDict_);
        Cmu_.readIfPresent(wallFunctionDict_);
        yPlusLam_ = yPlusLam(kappa_.value(), E_.value());

        k0_.readIfPresent(*this);
        epsilon0_.readIfPresent(*this);
        epsilonSmall_.readIfPresent(*this);

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
