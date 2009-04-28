/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

#include "faceAreaWeightModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(faceAreaWeightModel, 0);
defineRunTimeSelectionTable(faceAreaWeightModel, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

faceAreaWeightModel::faceAreaWeightModel
(
    const word& type,
    const dictionary& relaxationDict,
    const conformalVoronoiMesh& cvMesh
)
:
    dictionary(relaxationDict),
    cvMesh_(cvMesh),
    coeffDict_(subDict(type + "Coeffs"))
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<faceAreaWeightModel> faceAreaWeightModel::New
(
    const dictionary& relaxationDict,
    const conformalVoronoiMesh& cvMesh
)
{
    word faceAreaWeightModelTypeName
    (
        relaxationDict.lookup("faceAreaWeightModel")
    );

    Info<< nl << "Selecting faceAreaWeightModel "
        << faceAreaWeightModelTypeName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(faceAreaWeightModelTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "faceAreaWeightModel::New(const dictionary&, "
            "const conformalVoronoiMesh&)"
        )   << "Unknown faceAreaWeightModel type "
            << faceAreaWeightModelTypeName
            << endl << endl
            << "Valid faceAreaWeightModel types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<faceAreaWeightModel>(cstrIter()(relaxationDict, cvMesh));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

faceAreaWeightModel::~faceAreaWeightModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
