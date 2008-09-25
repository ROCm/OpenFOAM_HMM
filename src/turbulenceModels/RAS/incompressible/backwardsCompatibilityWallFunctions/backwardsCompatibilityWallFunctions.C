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

#include "backwardsCompatibilityWallFunctions.H"

#include "calculatedFvPatchField.H"
#include "nutWallFunctionFvPatchScalarField.H"
#include "epsilonWallFunctionFvPatchScalarField.H"
#include "kQRWallFunctionFvPatchField.H"
#include "omegaWallFunctionFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

wordList replaceWallBoundaryTypes
(
    const fvMesh& mesh,
    const wordList& oldTypeNames,
    const wordList& newTypeNames
)
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


tmp<volScalarField> autoCreateNut
(
    const word& fieldName,
    const fvMesh& mesh
)
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


tmp<volScalarField> autoCreateEpsilon
(
    const word& fieldName,
    const fvMesh& mesh
)
{
    return autoCreateWallFunctionField<scalar>
    (
        fieldName,
        mesh,
        RASModels::epsilonWallFunctionFvPatchScalarField::typeName
    );
}


tmp<volScalarField> autoCreateOmega
(
    const word& fieldName,
    const fvMesh& mesh
)
{
    return autoCreateWallFunctionField<scalar>
    (
        fieldName,
        mesh,
        RASModels::omegaWallFunctionFvPatchScalarField::typeName
    );
}


tmp<volScalarField> autoCreateK
(
    const word& fieldName,
    const fvMesh& mesh
)
{
    return autoCreateWallFunctionField<scalar>
    (
        fieldName,
        mesh,
        RASModels::kQRWallFunctionFvPatchField<scalar>::typeName
    );
}


tmp<volScalarField> autoCreateQ
(
    const word& fieldName,
    const fvMesh& mesh
)
{
    return autoCreateWallFunctionField<scalar>
    (
        fieldName,
        mesh,
        RASModels::kQRWallFunctionFvPatchField<scalar>::typeName
    );
}


tmp<volSymmTensorField> autoCreateR
(
    const word& fieldName,
    const fvMesh& mesh
)
{
    return autoCreateWallFunctionField<symmTensor>
    (
        fieldName,
        mesh,
        RASModels::kQRWallFunctionFvPatchField<symmTensor>::typeName
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //

