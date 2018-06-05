/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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

#include "proxyFvMesh.H"
#include "fvPatchField.H"

// * * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(proxyFvMesh, 0);
defineRunTimeSelectionTable(proxyFvMesh, time);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::proxyFvMesh::fvPatchFieldExists(const word& patchType)
{
    if
    (
        fvPatchField<scalar>::dictionaryConstructorTablePtr_->found(patchType)
     || fvPatchField<vector>::dictionaryConstructorTablePtr_->found(patchType)
     || fvPatchField<sphericalTensor>::
            dictionaryConstructorTablePtr_->found(patchType)
     || fvPatchField<symmTensor>::
            dictionaryConstructorTablePtr_->found(patchType)
     || fvPatchField<tensor>::dictionaryConstructorTablePtr_->found(patchType)
    )
    {
        return true;
    }

    return false;
}


Foam::proxyFvMesh::proxyFvMesh
(
    const IOobject& io,
    pointField&& points,
    faceList&& faces,
    labelList&& allOwner,
    labelList&& allNeighbour
)
:
    fvMesh
    (
        io,
        std::move(points),
        std::move(faces),
        std::move(allOwner),
        std::move(allNeighbour)
    )
{}


Foam::autoPtr<Foam::proxyFvMesh> Foam::proxyFvMesh::New
(
    const word& modelType,
    const Time& runTime
)
{
    Info<< "Selecting proxy mesh model " << modelType << endl;

    auto cstrIter = timeConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown proxy mesh type "
            << modelType << nl << nl
            << "Valid dumy meshes :" << endl
            << timeConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<proxyFvMesh>(cstrIter()(runTime));
}


// ************************************************************************* //
