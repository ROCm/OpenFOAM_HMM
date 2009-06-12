/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
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

#include "psiChemistryModel.H"

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::psiChemistryModel> Foam::psiChemistryModel::New
(
    const fvMesh& mesh
)
{
    word psiChemistryModelType;
    word thermoTypeName;
    word userSel;

    // Enclose the creation of the chemistrtyProperties to ensure it is
    // deleted before the chemistrtyProperties is created otherwise the
    // dictionary is entered in the database twice
    {
        IOdictionary chemistryPropertiesDict
        (
            IOobject
            (
                "chemistryProperties",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        chemistryPropertiesDict.lookup("psiChemistryModel") >> userSel;

        // construct chemistry model type name by inserting first template
        // argument
        label tempOpen = userSel.find('<');
        label tempClose = userSel.find('>');

        word className = userSel(0, tempOpen);
        thermoTypeName = userSel(tempOpen + 1, tempClose - tempOpen - 1);

        psiChemistryModelType =
            className + '<' + typeName + ',' + thermoTypeName + '>';
    }

    Info<< "Selecting psiChemistryModel " << userSel << endl;

    fvMeshConstructorTable::iterator cstrIter =
        fvMeshConstructorTablePtr_->find(psiChemistryModelType);

    if (cstrIter == fvMeshConstructorTablePtr_->end())
    {
        wordList models = fvMeshConstructorTablePtr_->toc();
        forAll(models, i)
        {
            models[i] = models[i].replace(typeName + ',', "");
        }

        FatalErrorIn("psiChemistryModelBase::New(const mesh&)")
            << "Unknown psiChemistryModel type " << userSel
            << nl << nl << "Valid psiChemistryModel types are:" << nl
            << models << nl << exit(FatalError);
    }

    return autoPtr<psiChemistryModel>
        (cstrIter()(mesh, typeName, thermoTypeName));
}


// ************************************************************************* //
