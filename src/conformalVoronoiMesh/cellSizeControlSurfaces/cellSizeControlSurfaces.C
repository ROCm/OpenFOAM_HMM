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

#include "cellSizeControlSurfaces.H"
#include "conformalVoronoiMesh.H"
#include "cellSizeFunction.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellSizeControlSurfaces::cellSizeControlSurfaces
(
    const conformalVoronoiMesh& cvMesh,
    const searchableSurfaces& allGeometry,
    const dictionary& motionControlDict
)
:
    cvMesh_(cvMesh),
    allGeometry_(allGeometry),
    surfaces_(),
    defaultCellSize_(readScalar(motionControlDict.lookup("defaultCellSize"))),
    defaultPriority_
    (
        motionControlDict.lookupOrDefault<label>("defaultPriority", 0)
    )
{
    const dictionary& surfacesDict
    (
        motionControlDict.subDict("cellSizeControlGeometry")
    );

    Info<< nl << "Reading cellSizeControlGeometry" << endl;

    surfaces_.setSize(surfacesDict.size());

    cellSizeFunctions_.setSize(surfacesDict.size());

    labelList priorities(surfacesDict.size());

    label surfI = 0;

    forAllConstIter(dictionary, surfacesDict, iter)
    {
        const dictionary& surfaceSubDict
        (
            surfacesDict.subDict(iter().keyword())
        );

        // If the "surface" keyword is not found in the dictionary, assume that
        // the name of the dictionary is the surface.  Distinction required to
        // allow the same surface to be used multiple times to supply multiple
        // cellSizeFunctions

        word surfaceName = surfaceSubDict.lookupOrDefault<word>
        (
            "surface",
            iter().keyword()
        );

        surfaces_[surfI] = allGeometry_.findSurfaceID(surfaceName);

        if (surfaces_[surfI] < 0)
        {
            FatalErrorIn
            (
                "Foam::cellSizeControlSurfaces::cellSizeControlSurfaces"
            )   << "No surface " << surfaceName << " found. "
                << "Valid geometry is " << nl << allGeometry_.names()
                << exit(FatalError);
        }

        const searchableSurface& surface = allGeometry_[surfaces_[surfI]];

        Info<< nl << "    " << iter().keyword() << nl
            << "    surface: " << surfaceName << endl;

        cellSizeFunctions_.set
        (
            surfI,
            cellSizeFunction::New
            (
                surfaceSubDict,
                cvMesh,
                surface
            )
        );

        priorities[surfI] = cellSizeFunctions_[surfI].priority();

        surfI++;
    }

    // Sort cellSizeFunctions_ and surfaces_ by priority.  Cut off any surfaces
    // where priority < defaultPriority_

    labelList sortedIndices;

    sortedOrder(priorities, sortedIndices);

    sortedIndices = invert(sortedIndices.size(), sortedIndices);

    // Reverse the sort order
    sortedIndices = (sortedIndices.size() - 1) - sortedIndices;

    inplaceReorder(sortedIndices, surfaces_);
    inplaceReorder(sortedIndices, priorities);
    cellSizeFunctions_.reorder(sortedIndices);

    forAll(priorities, surfI)
    {
        if (priorities[surfI] < defaultPriority_)
        {
            WarningIn("cellSizeControlSurfaces::cellSizeControlSurfaces")
                << "Priority of " << priorities[surfI]
                << " is less than defaultPriority " << defaultPriority_
                << ". All cellSizeFunctions with priorities lower than default "
                << "will be ignored."
                << endl;

            surfaces_.setSize(surfI);
            cellSizeFunctions_.setSize(surfI);

            break;
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellSizeControlSurfaces::~cellSizeControlSurfaces()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::cellSizeControlSurfaces::cellSize
(
    const point& pt,
    bool isSurfacePoint
) const
{
    scalar sizeAccumulator = 0;
    scalar numberOfFunctions = 0;

    label previousPriority =
        cellSizeFunctions_[cellSizeFunctions_.size() - 1].priority();

    forAll(cellSizeFunctions_, i)
    {
        const cellSizeFunction& cSF = cellSizeFunctions_[i];

        if (cSF.priority() < previousPriority && numberOfFunctions > 0)
        {
            return sizeAccumulator/numberOfFunctions;
        }

        scalar sizeI;

        if (cSF.cellSize(pt, sizeI, isSurfacePoint))
        {
            previousPriority = cSF.priority();

            sizeAccumulator += sizeI;
            numberOfFunctions++;
        }
    }

    if (previousPriority == defaultPriority_ || numberOfFunctions == 0)
    {
        sizeAccumulator += defaultCellSize_;
        numberOfFunctions++;
    }

    return sizeAccumulator/numberOfFunctions;
}


Foam::scalarField Foam::cellSizeControlSurfaces::cellSize
(
    const pointField& pts,
    const List<bool>& isSurfacePoint
) const
{
    if (pts.size() != isSurfacePoint.size())
    {
        FatalErrorIn
        (
            "Foam::cellSizeControlSurfaces::cellSizeControlSurfaces \
             ( \
                 const pointField& pts, \
                 const List<bool>& isSurfacePoint \
             ) \
             const"
        )   << "Size of pointField (" << pts.size()
            << ") and List<bool> (" << isSurfacePoint.size()
            << ") do not match." << nl
            << exit(FatalError);
    }

    scalarField cellSizes(pts.size());

    forAll(pts, i)
    {
        cellSizes[i] = cellSize(pts[i], isSurfacePoint[i]);
    }

    return cellSizes;
}



// ************************************************************************* //




