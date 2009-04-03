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

#include "initialPointsMethod.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(initialPointsMethod, 0);
defineRunTimeSelectionTable(initialPointsMethod, dictionary);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

initialPointsMethod::initialPointsMethod
(
    const word& type,
    const dictionary& initialPointsDict,
    const conformalVoronoiMesh& cvMesh
)
:
    dictionary(initialPointsDict),
    cvMesh_(cvMesh),
    detailsDict_(subDict(type + "Details"))
{}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<initialPointsMethod> initialPointsMethod::New
(
    const dictionary& initialPointsDict,
    const conformalVoronoiMesh& cvMesh
)
{
    word initialPointsMethodTypeName
    (
        initialPointsDict.lookup("initialPointsMethod")
    );

    Info<< nl << "Selecting initialPointsMethod "
        << initialPointsMethodTypeName << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(initialPointsMethodTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "initialPointsMethod::New(const volVectorField&, "
            "const surfaceScalarField&, transportModel&)"
        )   << "Unknown initialPointsMethod type "
            << initialPointsMethodTypeName
            << endl << endl
            << "Valid initialPointsMethod types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<initialPointsMethod>(cstrIter()(initialPointsDict, cvMesh));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

initialPointsMethod::~initialPointsMethod()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
