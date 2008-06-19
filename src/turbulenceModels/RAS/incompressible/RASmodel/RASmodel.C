/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

#include "RASmodel.H"
#include "wallFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(RASmodel, 0);
defineRunTimeSelectionTable(RASmodel, dictionary);

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void RASmodel::printCoeffs()
{
    if (printCoeffs_)
    {
        Info<< type() << "Coeffs" << RASmodelCoeffs_ << endl;;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

RASmodel::RASmodel
(
    const word& type,
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& lamTransportModel
)
:
    IOdictionary
    (
        IOobject
        (
            "turbulenceProperties",
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),

    runTime_(U.time()),
    mesh_(U.mesh()),

    U_(U),
    phi_(phi),
    transportModel_(lamTransportModel),

    turbulence_(lookup("turbulence")),
    printCoeffs_(lookupOrDefault<Switch>("printCoeffs", false)),
    RASmodelCoeffs_(subDict(type + "Coeffs")),

    kappa_
    (
        subDict("wallFunctionCoeffs").lookupOrAddDefault<scalar>
        (
            "kappa",
            0.4187
        )
    ),
    E_(subDict("wallFunctionCoeffs").lookupOrAddDefault<scalar>("E", 9.0)),
    yPlusLam_(yPlusLam(kappa_, E_)),

    k0_("k0", dimVelocity*dimVelocity, SMALL),
    epsilon0_("epsilon", k0_.dimensions()/dimTime, SMALL),
    epsilonSmall_("epsilonSmall", epsilon0_.dimensions(), SMALL),

    y_(mesh_)
{}


RASmodel::~RASmodel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

scalar RASmodel::yPlusLam(const scalar kappa, const scalar E)
{
    scalar ypl = 11.0;

    for (int i=0; i<10; i++)
    {
        ypl = log(E*ypl)/kappa;
    }

    return ypl;
}


tmp<scalarField> RASmodel::yPlus(const label patchNo) const
{
    const fvPatch& curPatch = mesh_.boundary()[patchNo];

    tmp<scalarField> tYp(new scalarField(curPatch.size()));
    scalarField& Yp = tYp();

    if (typeid(curPatch) == typeid(wallFvPatch))
    {
        scalar Cmu(readScalar(RASmodelCoeffs_.lookup("Cmu")));

        Yp = pow(Cmu, 0.25)*y_[patchNo]
            *sqrt(k()().boundaryField()[patchNo].patchInternalField())
            /nu().boundaryField()[patchNo];
    }
    else
    {
        WarningIn
        (
            "tmp<scalarField> RASmodel::yPlus(const label patchNo)"
        )   << "const : " << endl
            << "Patch " << patchNo << " is not a wall.  Returning blank field"
            << endl;

        Yp.setSize(0);
    }

    return tYp;
}


void RASmodel::correct()
{
    if (mesh_.changing())
    {
        y_.correct();
    }
}


bool RASmodel::read()
{
    if (regIOobject::read())
    {
        lookup("turbulence") >> turbulence_;
        RASmodelCoeffs_ = subDict(type() + "Coeffs");

        subDict("wallFunctionCoeffs").readIfPresent<scalar>("kappa", kappa_);
        subDict("wallFunctionCoeffs").readIfPresent<scalar>("E", E_);

        yPlusLam_ = yPlusLam(kappa_, E_);

        readIfPresent("k0", k0_);
        readIfPresent("epsilon0", epsilon0_);
        readIfPresent("epsilonSmall", epsilonSmall_);

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
