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

#include "cellSizeFunction.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

defineTypeNameAndDebug(cellSizeFunction, 0);
defineRunTimeSelectionTable(cellSizeFunction, dictionary);

scalar cellSizeFunction::snapToSurfaceTol_ = 1e-10;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

cellSizeFunction::cellSizeFunction
(
    const word& type,
    const dictionary& cellSizeFunctionDict,
    const conformalVoronoiMesh& cvMesh,
    const searchableSurface& surface
)
:
    dictionary(cellSizeFunctionDict),
    cvMesh_(cvMesh),
    surface_(surface),
    coeffsDict_(subDict(type + "Coeffs")),
    sideMode_(),
    priority_(readLabel(cellSizeFunctionDict.lookup("priority")))
{
    word mode = cellSizeFunctionDict.lookup("mode");

    if (surface_.hasVolumeType())
    {
        if (mode == "inside")
        {
            sideMode_ = INSIDE;
        }
        else if (mode == "outside")
        {
            sideMode_ = OUTSIDE;
        }
        else if (mode == "bothSides")
        {
            sideMode_ = BOTHSIDES;
        }
        else
        {
            FatalErrorIn("cellSizeFunction::cellSizeFunction")
            << "Unknown mode, expected: inside, outside or bothSides" << nl
                << exit(FatalError);
        }
    }
    else if (mode != BOTHSIDES)
    {
        WarningIn("cellSizeFunction::cellSizeFunction")
            << "surface does not support volumeType, defaulting mode to "
            << "bothSides."
            << endl;

        sideMode_ = BOTHSIDES;
    }
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<cellSizeFunction> cellSizeFunction::New
(
    const dictionary& cellSizeFunctionDict,
    const conformalVoronoiMesh& cvMesh,
    const searchableSurface& surface
)
{
    word cellSizeFunctionTypeName
    (
        cellSizeFunctionDict.lookup("cellSizeFunction")
    );

    Info<< "    Selecting cellSizeFunction " << cellSizeFunctionTypeName
        << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(cellSizeFunctionTypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "cellSizeFunction::New(dictionary&, "
            "const conformalVoronoiMesh&, const searchableSurface&)"
        )   << "Unknown cellSizeFunction type "
            << cellSizeFunctionTypeName
            << endl << endl
            << "Valid cellSizeFunction types are :" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<cellSizeFunction>
    (
        cstrIter()(cellSizeFunctionDict, cvMesh, surface)
    );
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

cellSizeFunction::~cellSizeFunction()
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
