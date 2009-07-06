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
        motionControlDict.lookupOrDefault<scalar>("defaultPriority", 0)
    )
{
    const dictionary& surfacesDict
    (
        motionControlDict.subDict("cellSizeControlGeometry")
    );

    Info<< nl << "Reading cellSizeControlGeometry." << endl;

    surfaces_.setSize(surfacesDict.size());

    label surfI = 0;

    forAllConstIter(dictionary, surfacesDict, iter)
    {
        word surfaceName = iter().keyword();

        surfaces_[surfI] = allGeometry_.findSurfaceID(surfaceName);

        if (surfaces_[surfI] < 0)
        {
            FatalErrorIn
            (
                "Foam::cellSizeControlSurfaces::cellSizeControlSurfaces"
            )   << "No surface " << iter().keyword() << " found. "
                << "Valid geometry is " << nl << allGeometry_.names()
                << exit(FatalError);
        }

        const dictionary& surfaceSubDict(surfacesDict.subDict(surfaceName));

        // Code here ->

        Info<< nl << surfaceName << surfaceSubDict << endl;

        // <- Code here

        surfI++;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellSizeControlSurfaces::~cellSizeControlSurfaces()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::scalar Foam::cellSizeControlSurfaces::targetCellSize
(
    const point& pt
) const
{
    Info<< "TARGETCELLSIZE NOT IMPLEMENTED." << endl;

    return 1.0;
}


// ************************************************************************* //




