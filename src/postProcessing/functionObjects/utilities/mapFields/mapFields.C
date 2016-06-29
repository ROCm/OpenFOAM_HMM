/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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

#include "mapFields.H"
#include "meshToMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mapFieldsFO, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::mapFieldsFO::createInterpolation(const dictionary& dict)
{
    const fvMesh& meshTarget = static_cast<const fvMesh&>(obr_);
    const word mapRegionName(dict.lookup("mapRegion"));

    if (log_)
    {
        Info<< name_ << ':' << nl
            << "    Reading mesh " << mapRegionName << endl;
    }

    mapRegionPtr_.reset
    (
        new fvMesh
        (
            IOobject
            (
                mapRegionName,
                meshTarget.time().constant(),
                meshTarget.time()
            )
        )
    );
    const fvMesh& mapRegion = mapRegionPtr_();
    word mapMethodName(dict.lookup("mapMethod"));
    if (!meshToMesh::interpolationMethodNames_.found(mapMethodName))
    {
        FatalErrorInFunction
            << type() << " " << name_ << ": unknown map method "
            << mapMethodName << nl
            << "Available methods include: "
            << meshToMesh::interpolationMethodNames_.sortedToc()
            << exit(FatalError);
    }

    meshToMesh::interpolationMethod mapMethod
    (
        meshToMesh::interpolationMethodNames_[mapMethodName]
    );

    // Lookup corresponding AMI method
    word patchMapMethodName =
        AMIPatchToPatchInterpolation::interpolationMethodToWord
        (
            meshToMesh::interpolationMethodAMI(mapMethod)
        );

    // Optionally override
    if (dict.readIfPresent("patchMapMethod", patchMapMethodName))
    {
        if (log_)
            Info<< "    Patch mapping method: " << patchMapMethodName << endl;
    }

    bool consistent = readBool(dict.lookup("consistent"));

    if (log_) Info<< "    Creating mesh to mesh interpolation" << endl;

    if (consistent)
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
        HashTable<word> patchMap(dict.lookup("patchMap"));
        wordList cuttingPatches(dict.lookup("cuttingPatches"));

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

Foam::mapFieldsFO::mapFieldsFO
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    active_(true),
    log_(true),
    mapRegionPtr_(),
    interpPtr_(),
    fieldNames_()
{
    if (!isA<fvMesh>(obr_))
    {
        active_ = false;
        WarningInFunction
            << "No fvMesh available, deactivating " << name_ << nl
            << endl;
    }

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mapFieldsFO::~mapFieldsFO()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mapFieldsFO::read(const dictionary& dict)
{
    if (active_)
    {
        log_.readIfPresent("log", dict);
        dict.lookup("fields") >> fieldNames_;

        createInterpolation(dict);
    }
}


void Foam::mapFieldsFO::execute()
{}


void Foam::mapFieldsFO::end()
{}


void Foam::mapFieldsFO::timeSet()
{}


void Foam::mapFieldsFO::write()
{
    if (active_)
    {
        if (log_) Info
            << type() << " " << name_ << " output:" << nl;

        bool ok = false;

        ok = writeFieldType<scalar>() || ok;
        ok = writeFieldType<vector>() || ok;
        ok = writeFieldType<sphericalTensor>() || ok;
        ok = writeFieldType<symmTensor>() || ok;
        ok = writeFieldType<tensor>() || ok;

        if (log_)
        {
            if (!ok)
            {
                Info<< "    none" << nl;
            }

            Info<< endl;
        }
    }
}


// ************************************************************************* //
