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
#include "kQRWallFunctionFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
RASModel::autoCreateWallFunctionField
(
    const word& fieldName,
    const fvMesh& mesh,
    const word& wallFunctionName
) const
{
    IOobject nutHeader
    (
        "nut",
        mesh.time().timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    if (nutHeader.headerOk())
    {
        return tmp<fieldType>
        (
            new fieldType
            (
                IOobject
                (
                    fieldName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh
            )
        );
    }
    else
    {
        Info<< "--> Upgrading " << fieldName << " to employ run-time "
            << "selectable wall functions" << endl;

        // Read existing epsilon field
        tmp<fieldType> fieldOrig
        (
            new fieldType
            (
                IOobject
                (
                    fieldName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh
            )
        );

        wordList fieldBoundaryTypes = replaceWallBoundaryTypes
        (
            mesh,
            fieldOrig().boundaryField().types(),
            wordList
            (
                fieldOrig().boundaryField().types().size(),
                wallFunctionName
            )
        );

        tmp<fieldType> fieldNew
        (
            new fieldType
            (
                IOobject
                (
                    fieldName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensioned<Type>
                (
                    "zero",
                    fieldOrig().dimensions(),
                    pTraits<Type>::zero
                ),
                fieldBoundaryTypes
            )
        );

        fieldNew() == fieldOrig();

        Info<< "    Writing backup of original " << fieldName << " to "
            << fieldName << ".old" << endl;
        fieldOrig().rename(fieldName + ".old");
        fieldOrig().write();

        Info<< "    Writing updated " << fieldName << endl;
        fieldNew().write();

        return fieldNew;
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> > RASModel::autoCreateKQR
(
    const word& fieldName,
    const fvMesh& mesh
) const
{
    return autoCreateWallFunctionField<Type>
    (
        fieldName,
        mesh,
        RASModels::kQRWallFunctionFvPatchField<Type>::typeName
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
