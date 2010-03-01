/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

Description
    Set normal consistent with respect to a user provided 'outside' point.
    If -inside the point is considered inside.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "orientedSurface.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();
    argList::validArgs.append("Foam surface file");
    argList::validArgs.append("visiblePoint");
    argList::validArgs.append("output file");
    argList::addBoolOption("inside");
    argList args(argc, argv);

    const fileName surfFileName = args[1];
    const point visiblePoint    = args.argRead<point>(2);
    const fileName outFileName  = args[3];

    const bool orientInside = args.optionFound("inside");

    Info<< "Reading surface from " << surfFileName << endl;
    Info<< "Visible point " << visiblePoint << endl;

    if (orientInside)
    {
        Info<< "Orienting surface such that visiblePoint " << visiblePoint
            << " is inside" << endl;
    }
    else
    {
        Info<< "Orienting surface such that visiblePoint " << visiblePoint
            << " is outside" << endl;
    }

    Info<< "Writing surface to " << outFileName << endl;


    // Load surface
    triSurface surf(surfFileName);

    //orientedSurface normalSurf(surf, visiblePoint, !orientInside);
    bool anyFlipped = orientedSurface::orient
    (
        surf,
        visiblePoint,
       !orientInside
    );

    if (anyFlipped)
    {
        Info<< "Flipped orientation of (part of) surface." << endl;
    }
    else
    {
        Info<< "Did not flip orientation of any triangle of surface." << endl;
    }

    Info<< "Writing new surface to " << outFileName << endl;

    surf.write(outFileName);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
