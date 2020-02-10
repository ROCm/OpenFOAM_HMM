/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2020 OpenCFD Ltd.
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

InNamespace
    Foam

Description
    Helper routines for reading a field or fields,
    for foamToEnsight

SourceFiles
    readFields.C

\*---------------------------------------------------------------------------*/

#ifndef readFields_H
#define readFields_H

#include "instantList.H"
#include "IOobjectList.H"
#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Get the field or return nullptr
template<class GeoField>
tmp<GeoField> getField
(
    const IOobject* io,
    const typename GeoField::Mesh& mesh
)
{
    if (io)
    {
        return tmp<GeoField>::New(*io, mesh);
    }

    return nullptr;
}


//- Get the named field from the objects, or return nullptr.
template<class GeoField>
tmp<GeoField> getField
(
    const typename GeoField::Mesh& mesh,
    const IOobjectList& objects,
    const word& fieldName
)
{
    // Can do something with syncPar on failure ...

    return getField<GeoField>(objects.findObject(fieldName), mesh);
}


//- Convert an internal field to zero-gradient volume field
template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
makeZeroGradientField
(
    const tmp
    <
        typename GeometricField<Type, fvPatchField, volMesh>::Internal
    >& tdf
)
{
    if (tdf.valid())
    {
        auto& df = tdf.ref();

        auto tfield = GeometricField<Type, fvPatchField, volMesh>::New
        (
            df.name(),
            df.mesh(),
            df.dimensions(),
            std::move(df.field()),
            zeroGradientFvPatchScalarField::typeName
        );

        tfield.ref().oriented() = df.oriented();
        tfield.ref().correctBoundaryConditions();

        tdf.clear();

        return tfield;
    }

    tdf.clear();

    return nullptr;
}


//- Convert a volume field to zero-gradient volume field
template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
makeZeroGradientField
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tdf
)
{
    if (tdf.valid())
    {
        auto& df = tdf.ref();

        auto tfield = GeometricField<Type, fvPatchField, volMesh>::New
        (
            df.name(),
            df.mesh(),
            df.dimensions(),
            std::move(df.primitiveFieldRef(false)),  // No update accessTime
            zeroGradientFvPatchScalarField::typeName
        );

        tfield.ref().oriented() = df.oriented();
        tfield.ref().correctBoundaryConditions();

        tdf.clear();

        return tfield;
    }

    tdf.clear();

    return nullptr;
}


//- Check if fields are good to use (available at all times)
//  ignore special fields (_0 fields),
//  ignore fields that are not available for all time-steps
label checkData
(
    const fvMesh& mesh,
    const instantList& timeDirs,
    wordList& objectNames
);


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
