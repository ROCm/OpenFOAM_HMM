/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2011 OpenCFD Ltd.
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

#include "cellSizeControlSurfaces.H"
#include "conformalVoronoiMesh.H"
#include "cellSizeFunction.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(cellSizeControlSurfaces, 0);

}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

bool Foam::cellSizeControlSurfaces::evalCellSizeFunctions
(
    const point& pt,
    scalar& minSize,
    bool isSurfacePoint
) const
{
    bool anyFunctionFound = false;

    // // Regions requesting with the same priority take the average

    // scalar sizeAccumulator = 0;
    // scalar numberOfFunctions = 0;

    // label previousPriority = defaultPriority_;

    // if (cellSizeFunctions_.size())
    // {
    //     previousPriority =
    //         cellSizeFunctions_[cellSizeFunctions_.size() - 1].priority();

    //     forAll(cellSizeFunctions_, i)
    //     {
    //         const cellSizeFunction& cSF = cellSizeFunctions_[i];

    //         if (cSF.priority() < previousPriority && numberOfFunctions > 0)
    //         {
    //             return sizeAccumulator/numberOfFunctions;
    //         }

    //         scalar sizeI;

    //         if (cSF.cellSize(pt, sizeI, isSurfacePoint))
    //         {
    //             anyFunctionFound = true;

    //             previousPriority = cSF.priority();

    //             sizeAccumulator += sizeI;
    //             numberOfFunctions++;
    //         }
    //     }
    // }

    // if (previousPriority == defaultPriority_ || numberOfFunctions == 0)
    // {
    //     sizeAccumulator += defaultCellSize_;
    //     numberOfFunctions++;
    // }

    // minSize = sizeAccumulator/numberOfFunctions;

    // return anyFunctionFound;

    // Regions requesting with the same priority take the smallest

    if (cellSizeFunctions_.size())
    {
        // Initialise to the last (lowest) priority
        label previousPriority = cellSizeFunctions_.last().priority();

        forAll(cellSizeFunctions_, i)
        {
            const cellSizeFunction& cSF = cellSizeFunctions_[i];

            if (debug)
            {
                Info<< "size function "
                    << allGeometry_.names()[surfaces_[i]]
                    << " priority " << cSF.priority()
                    << endl;
            }

            if (cSF.priority() < previousPriority)
            {
                return minSize;
            }

            scalar sizeI;

            if (cSF.cellSize(pt, sizeI, isSurfacePoint))
            {
                anyFunctionFound = true;

                if (cSF.priority() == previousPriority)
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
                    Info<< "sizeI " << sizeI << " minSize " << minSize << endl;
                }

                previousPriority = cSF.priority();
            }
        }
    }

    return anyFunctionFound;
}


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
    cellSizeFunctions_(),
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
    scalar size = defaultCellSize_;

    bool anyFunctionFound = evalCellSizeFunctions(pt, size, isSurfacePoint);

    if (!anyFunctionFound)
    {
        // Check if the point in question was actually inside the domain, if
        // not, then it may be falling back to an inappropriate default size.

        if (cvMesh_.geometryToConformTo().outside(pt))
        {
            pointIndexHit surfHit;
            label hitSurface;

            cvMesh_.geometryToConformTo().findSurfaceNearest
            (
                pt,
                cvMesh_.geometryToConformTo().spanMagSqr(),
                surfHit,
                hitSurface
            );

            if (!surfHit.hit())
            {
                FatalErrorIn
                (
                    "Foam::scalar Foam::cellSizeControlSurfaces::cellSize"
                    "("
                        "const point& pt, "
                        "bool isSurfacePoint"
                    ") const"
                )
                    << "Point " << pt << "was not within "
                    << cvMesh_.geometryToConformTo().spanMag()
                    << " of the surface." << nl
                    << "findSurfaceNearest did not find a hit across the span "
                    << "of the surfaces."
                    << nl << exit(FatalError) << endl;
            }
            else
            {
                // Evaluating the cell size at the nearest surface
                evalCellSizeFunctions(surfHit.hitPoint(), size, true);
            }
        }
    }

    return size;
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
