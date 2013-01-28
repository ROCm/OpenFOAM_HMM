/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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

#include "cellSizeAndAlignmentControls.H"
#include "searchableSurfaceControl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(cellSizeAndAlignmentControls, 0);
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

bool Foam::cellSizeAndAlignmentControls::evalCellSizeFunctions
(
    const point& pt,
    scalar& minSize
) const
{
    bool anyFunctionFound = false;

    // Regions requesting with the same priority take the smallest

    if (controlFunctions_.size())
    {
        // Maintain priority of current hit. Initialise so it always goes
        // through at least once.
        label previousPriority = -1;

        forAll(controlFunctions_, i)
        {
            const cellSizeAndAlignmentControl& cSF = controlFunctions_[i];

            if (isA<searchableSurfaceControl>(cSF))
            {
                const searchableSurfaceControl& sSC =
                    refCast<const searchableSurfaceControl>(cSF);

                if (debug)
                {
                    Info<< "size function "
                        << sSC.name()
                        << " priority " << sSC.priority()
                        << endl;
                }

                if (sSC.priority() < previousPriority)
                {
                    return minSize;
                }

                scalar sizeI;

                if (sSC.sizeFunction().cellSize(pt, sizeI))
                {
                    anyFunctionFound = true;

                    if (sSC.priority() == previousPriority)
                    {
                        if (sizeI < minSize)
                        {
                            minSize = sizeI;
                        }
                    }
                    else
                    {
                        minSize = sizeI;
                    }

                    if (debug)
                    {
                        Info<< "sizeI " << sizeI
                            <<" minSize " << minSize << endl;
                    }

                    previousPriority = sSC.priority();
                }
            }
        }
    }

    return anyFunctionFound;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellSizeAndAlignmentControls::cellSizeAndAlignmentControls
(
    const Time& runTime,
    const dictionary& shapeControlDict,
    const conformationSurfaces& geometryToConformTo,
    const scalar defaultCellSize
)
:
    shapeControlDict_(shapeControlDict),
    geometryToConformTo_(geometryToConformTo),
    controlFunctions_(shapeControlDict_.size()),
    defaultCellSize_(defaultCellSize)
{
    label functionI = 0;

    forAllConstIter(dictionary, shapeControlDict_, iter)
    {
        word shapeControlEntryName = iter().keyword();

        const dictionary& controlFunctionDict
        (
            shapeControlDict_.subDict(shapeControlEntryName)
        );

        Info<< nl << "Shape Control : " << shapeControlEntryName << endl;

        controlFunctions_.set
        (
            functionI,
            cellSizeAndAlignmentControl::New
            (
                runTime,
                shapeControlEntryName,
                controlFunctionDict,
                geometryToConformTo
            )
        );

        functionI++;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellSizeAndAlignmentControls::~cellSizeAndAlignmentControls()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::cellSizeAndAlignmentControls::cellSize
(
    const point& pt
) const
{
    scalar size = defaultCellSize_;

    evalCellSizeFunctions(pt, size);

    return size;
}


// ************************************************************************* //
