/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "mapFields.H"
#include "meshToMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(mapFields, 0);
    addToRunTimeSelectionTable(functionObject, mapFields, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::mapFields::createInterpolation
(
    const dictionary& dict
)
{
    const fvMesh& meshTarget = mesh_;
    const word mapRegionName(dict.get<word>("mapRegion"));

    Info<< name() << ':' << nl
        << "    Reading mesh " << mapRegionName << endl;

    mapRegionPtr_.reset
    (
        new fvMesh
        (
            IOobject
            (
                mapRegionName,
                meshTarget.time().constant(),
                meshTarget.time(),
                IOobject::MUST_READ
           )
        )
    );

    const fvMesh& mapRegion = mapRegionPtr_();
    const word mapMethodName(dict.get<word>("mapMethod"));
    if (!meshToMesh::interpolationMethodNames_.found(mapMethodName))
    {
        FatalErrorInFunction
            << type() << " " << name() << ": unknown map method "
            << mapMethodName << nl
            << "Available methods include: "
            << meshToMesh::interpolationMethodNames_
            << exit(FatalError);
    }

    meshToMesh::interpolationMethod mapMethod
    (
        meshToMesh::interpolationMethodNames_[mapMethodName]
    );

    // Lookup corresponding AMI method
    word patchMapMethodName = meshToMesh::interpolationMethodAMI(mapMethod);

    // Optionally override
    if (dict.readIfPresent("patchMapMethod", patchMapMethodName))
    {
        Info<< "    Patch mapping method: " << patchMapMethodName << endl;
    }

    Info<< "    Creating mesh to mesh interpolation" << endl;

    if (dict.get<bool>("consistent"))
    {
        interpPtr_.reset
        (
            new meshToMesh
            (
                mapRegion,
                meshTarget,
                mapMethodName,
                patchMapMethodName
            )
        );
    }
    else
    {
        HashTable<word> patchMap;
        wordList cuttingPatches;

        dict.readEntry("patchMap", patchMap);
        dict.readEntry("cuttingPatches", cuttingPatches);

        interpPtr_.reset
        (
            new meshToMesh
            (
                mapRegion,
                meshTarget,
                mapMethodName,
                patchMapMethodName,
                patchMap,
                cuttingPatches
            )
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::mapFields::mapFields
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    mapRegionPtr_(),
    interpPtr_(),
    fieldNames_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::mapFields::read(const dictionary& dict)
{
    if (fvMeshFunctionObject::read(dict))
    {
        dict.readEntry("fields", fieldNames_);
        createInterpolation(dict);

        return true;
    }

    return false;
}


bool Foam::functionObjects::mapFields::execute()
{
    Log << type() << " " << name() << " execute:" << nl;

    bool ok = false;

    ok = mapFieldType<scalar>() || ok;
    ok = mapFieldType<vector>() || ok;
    ok = mapFieldType<sphericalTensor>() || ok;
    ok = mapFieldType<symmTensor>() || ok;
    ok = mapFieldType<tensor>() || ok;

    if (log)
    {
        if (!ok)
        {
            Info<< "    none" << nl;
        }

        Info<< endl;
    }
    return true;
}


bool Foam::functionObjects::mapFields::write()
{
    Log << type() << " " << name() << " write:" << nl;

    bool ok = false;

    ok = writeFieldType<scalar>() || ok;
    ok = writeFieldType<vector>() || ok;
    ok = writeFieldType<sphericalTensor>() || ok;
    ok = writeFieldType<symmTensor>() || ok;
    ok = writeFieldType<tensor>() || ok;

    if (log)
    {
        if (!ok)
        {
            Info<< "    none" << nl;
        }

        Info<< endl;
    }

    return true;
}


// ************************************************************************* //
