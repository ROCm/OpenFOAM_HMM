/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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

#include "multiphaseSystem.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::multiphaseSystem> Foam::multiphaseSystem::New
(
    const fvMesh& mesh
)
{
    const word multiphaseSystemType
    (
        IOdictionary
        (
            IOobject
            (
                phasePropertiesName,
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        ).get<word>("type")
    );

    Info<< "Selecting multiphaseSystem " << multiphaseSystemType << endl;

    const auto cstrIter =
        dictionaryConstructorTablePtr_->cfind(multiphaseSystemType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown multiphaseSystemType type "
            << multiphaseSystemType << endl
            << "Valid multiphaseSystem types are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<multiphaseSystem> (cstrIter()(mesh));
}


// ************************************************************************* //
