/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
     \\/     M anipulation  | Copyright (C) 2015 IH-Cantabria
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

#include "waveModel.H"
#include "fvMesh.H"

Foam::autoPtr<Foam::waveModel> Foam::waveModel::New
(
    const word& dictName,
    const fvMesh& mesh,
    const polyPatch& patch
)
{
    IOdictionary waveDict
    (
        IOobject
        (
            dictName,
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false                   // Not registering
        )
    );

    word modelType = "none";
    dictionary patchDict;
    if (waveDict.found(patch.name()))
    {
        patchDict = waveDict.subDict(patch.name());
        patchDict.lookup("waveModel") >> modelType;
    }
    else
    {
        FatalIOErrorInFunction(waveDict)
            << "Dictionary entry for patch " << patch.name() << " not found"
            << exit(FatalIOError);
    }

    Info<< "Selecting waveModel " << modelType << endl;

    auto cstrIter = patchConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.found())
    {
        FatalIOErrorInFunction(waveDict)
            << "Unknown waveModel type "
            << modelType << nl << nl
            << "Valid waveModel types :" << nl
            << patchConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    return autoPtr<waveModel>(cstrIter()(patchDict, mesh, patch));
}


Foam::tmp<Foam::waveModel> Foam::waveModel::lookupOrCreate
(
    const polyPatch& patch,
    const fvMesh& mesh,
    const word& waveDictName
)
{
    const word modelName = waveModel::modelName(patch.name());

    if (!mesh.foundObject<waveModel>(modelName))
    {
        autoPtr<waveModel> model(waveModel::New(waveDictName, mesh, patch));
        waveModel* waveModelPtr = model.ptr();
        waveModelPtr->store();
        waveModelPtr->info(Info);
    }

    return mesh.lookupObject<waveModel>(modelName);
}


// ************************************************************************* //
