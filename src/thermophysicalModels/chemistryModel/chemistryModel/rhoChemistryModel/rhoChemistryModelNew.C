/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
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

#include "rhoChemistryModel.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::rhoChemistryModel>
Foam::rhoChemistryModel::New
(
    const fvMesh& mesh
)
{
    // get model name, but do not register the dictionary
    // otherwise it is registered in the database twice
    const word userModel
    (
        IOdictionary
        (
            IOobject
            (
                "chemistryProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        ).lookup("rhoChemistryModel")
    );

    // construct chemistry model type name by inserting first template argument
    const label tempOpen = userModel.find('<');
    const label tempClose = userModel.find('>');

    const word className = userModel(0, tempOpen);
    const word thermoTypeName =
        userModel(tempOpen + 1, tempClose - tempOpen - 1);

    const word modelType =
        className + '<' + typeName + ',' + thermoTypeName + '>';

    if (debug)
    {
        Info<< "Selecting rhoChemistryModel " << modelType << endl;
    }
    else
    {
        Info<< "Selecting rhoChemistryModel " << userModel << endl;
    }

    fvMeshConstructorTable::iterator cstrIter =
        fvMeshConstructorTablePtr_->find(modelType);

    if (cstrIter == fvMeshConstructorTablePtr_->end())
    {
        if (debug)
        {
            FatalErrorIn("rhoChemistryModelBase::New(const mesh&)")
                << "Unknown rhoChemistryModel type "
                << modelType << nl << nl
                << "Valid rhoChemistryModel types are:" << nl
                << fvMeshConstructorTablePtr_->sortedToc() << nl
                << exit(FatalError);
        }
        else
        {
            wordList models = fvMeshConstructorTablePtr_->sortedToc();
            forAll(models, i)
            {
                models[i] = models[i].replace(typeName + ',', "");
            }

            FatalErrorIn("rhoChemistryModelBase::New(const mesh&)")
                << "Unknown rhoChemistryModel type "
                << userModel << nl << nl
                << "Valid rhoChemistryModel types are:" << nl
                << models << nl
                << exit(FatalError);
        }
    }

    return autoPtr<rhoChemistryModel>
        (cstrIter()(mesh, typeName, thermoTypeName));
}


// ************************************************************************* //
