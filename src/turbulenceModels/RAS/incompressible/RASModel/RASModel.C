/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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

// Headers reqd for back-porting wall function boundary conditions
#include "calculatedFvPatchField.H"
#include "nutWallFunctionFvPatchScalarField.H"
#include "epsilonWallFunctionFvPatchScalarField.H"
#include "omegaWallFunctionFvPatchScalarField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(RASModel, 0);
defineRunTimeSelectionTable(RASModel, dictionary);

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void RASModel::printCoeffs()
{
    if (printCoeffs_)
    {
        Info<< type() << "Coeffs" << coeffDict_ << endl;;
    }
}


wordList RASModel::replaceWallBoundaryTypes
(
    const fvMesh& mesh,
    const wordList& oldTypeNames,
    const wordList& newTypeNames
) const
{
    const fvBoundaryMesh& bm = mesh.boundary();

    wordList boundaryTypes(bm.size());

    forAll(bm, patchI)
    {
        if (isType<wallFvPatch>(bm[patchI]))
        {
            boundaryTypes[patchI] = newTypeNames[patchI];
        }
        else
        {
            boundaryTypes[patchI] = oldTypeNames[patchI];
        }
    }

    return boundaryTypes;
}


tmp<volScalarField> RASModel::autoCreateNut
(
    const word& fieldName,
    const fvMesh& mesh
) const
{
    IOobject nutHeader
    (
        fieldName,
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    if (nutHeader.headerOk())
    {
        return tmp<volScalarField>(new volScalarField(nutHeader, mesh));
    }
    else
    {
        Info<< "--> Upgrading " << fieldName << " to employ run-time "
            << "selectable wall functions" << endl;

        wordList nutBoundaryTypes = replaceWallBoundaryTypes
        (
            mesh,
            wordList
            (
                mesh.boundary().size(),
                calculatedFvPatchField<scalar>::typeName
            ),
            wordList
            (
                mesh.boundary().size(),
                RASModels::nutWallFunctionFvPatchScalarField::typeName
            )
        );

        tmp<volScalarField> nut
        (
            new volScalarField
            (
                IOobject
                (
                    fieldName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("zero", dimArea/dimTime, 0.0),
                nutBoundaryTypes
            )
        );

        Info<< "    Writing updated " << fieldName << endl;
        nut().write();

        return nut;
    }
}


tmp<volScalarField> RASModel::autoCreateEpsilon
(
    const word& fieldName,
    const fvMesh& mesh
) const
{
    return autoCreateWallFunctionField<scalar>
    (
        fieldName,
        mesh,
        RASModels::epsilonWallFunctionFvPatchScalarField::typeName
    );
}


tmp<volScalarField> RASModel::autoCreateOmega
(
    const word& fieldName,
    const fvMesh& mesh
) const
{
    return autoCreateWallFunctionField<scalar>
    (
        fieldName,
        mesh,
        RASModels::omegaWallFunctionFvPatchScalarField::typeName
    );
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

    runTime_(U.time()),
    mesh_(U.mesh()),

    U_(U),
    phi_(phi),
    transportModel_(lamTransportModel),

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

    if (typeid(curPatch) == typeid(wallFvPatch))
    {
        Yp = pow(Cmu_.value(), 0.25)*y_[patchNo]
            *sqrt(k()().boundaryField()[patchNo].patchInternalField())
            /nu().boundaryField()[patchNo];
    }
    else
    {
        WarningIn
        (
            "tmp<scalarField> RASModel::yPlus(const label patchNo)"
        )   << "const : " << nl
            << "Patch " << patchNo << " is not a wall.  Returning zero field"
            << nl << endl;

        Yp.setSize(0);
    }

    return tYp;
}


void RASModel::correct()
{
    if (mesh_.changing())
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
