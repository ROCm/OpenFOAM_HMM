/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2011 OpenCFD Ltd.
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

Application
    surfaceMeshInfo

Description
    Miscellaneous information about surface meshes

Usage
    - surfaceMeshInfo surfaceFile [OPTION]

    \param -areas \n
    Report area for each face.

    \param -scale \<scale\> \n
    Specify a scaling factor when reading files.

    \param -xml \n
    Write output in XML format.

Note
    The filename extensions are used to determine the file format type.

    The XML-like output can be useful for extraction with other tools,
    but either output format can be easily extracted with a simple sed
    command:
    \verbatim
        surfaceMeshInfo surfaceFile -areas | \
            sed -ne '/areas/,/:/{ /:/!p }'

        surfaceMeshInfo surfaceFile -areas -xml | \
            sed -ne '/<areas/,/</{ /</!p }'
    \endverbatim

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "Time.H"

#include "UnsortedMeshedSurfaces.H"

using namespace Foam;


//
// duplicates code from primitiveMeshFaceCentresAndAreas.C
// should be refactored
//
Foam::vector faceArea(const face& fs, const pointField& p)
{
    const labelList& f = fs;
    label nPoints = f.size();

    point fCentre;
    vector fArea;

    // If the face is a triangle, do a direct calculation for efficiency
    // and to avoid round-off error-related problems
    if (nPoints == 3)
    {
        fCentre = (1.0/3.0)*(p[f[0]] + p[f[1]] + p[f[2]]);
        fArea = 0.5*((p[f[1]] - p[f[0]])^(p[f[2]] - p[f[0]]));
    }
    else
    {
        fCentre = p[f[0]];
        for (label pi = 1; pi < nPoints; ++pi)
        {
            fCentre += p[f[pi]];
        }

        fCentre /= nPoints;

        vector sumN = vector::zero;
        scalar sumA = 0.0;
        vector sumAc = vector::zero;

        for (label pi = 0; pi < nPoints; ++pi)
        {
            const point& nextPoint = p[f[(pi + 1) % nPoints]];

            const vector c = p[f[pi]] + nextPoint + fCentre;
            const vector n = (nextPoint - p[f[pi]])^(fCentre - p[f[pi]]);
            const scalar a = mag(n);

            sumN += n;
            sumA += a;
            sumAc += a*c;
        }

        fArea = 0.5*sumN;
    }

    return fArea;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "information about surface meshes"
    );

    argList::noBanner();
    argList::noParallel();
    argList::validArgs.append("surfaceFile");

    argList::addOption
    (
        "scale",
        "factor",
        "geometry scaling factor - default is 1"
    );
    argList::addBoolOption
    (
        "areas",
        "display area of each face"
    );
    argList::addBoolOption
    (
        "xml",
        "write output in XML format"
    );

    argList args(argc, argv);
    Time runTime(args.rootPath(), args.caseName());

    const fileName importName = args[1];

    // check that reading is supported
    if (!UnsortedMeshedSurface<face>::canRead(importName, true))
    {
        return 1;
    }

    const bool writeXML = args.optionFound("xml");
    const bool writeAreas = args.optionFound("areas");


    // use UnsortedMeshedSurface, not MeshedSurface to maintain ordering
    UnsortedMeshedSurface<face> surf(importName);

    scalar scaling = 0;
    if (args.optionReadIfPresent("scale", scaling) && scaling > 0)
    {
        Info<< " -scale " << scaling << endl;
        surf.scalePoints(scaling);
    }

    scalar areaTotal = 0;

    if (writeXML)
    {
        Info<<"<?xml version='1.0' encoding='utf-8'?>" << nl
            <<"<surfaceMeshInfo>" << nl
            << "<npoints>" << surf.nPoints() << "</npoints>" << nl
            << "<nfaces>"  << surf.size() << "</nfaces>" << nl;

        if (writeAreas)
        {
            Info<<"<areas size='" << surf.size() << "'>" << nl;
        }
    }
    else
    {
        Info<< "nPoints   : " << surf.nPoints() << nl
            << "nFaces    : " << surf.size() << nl;

        if (writeAreas)
        {
            Info<< "areas     : " << nl;
        }
    }

    forAll(surf, faceI)
    {
        scalar fArea = mag(faceArea(surf[faceI], surf.points()));
        areaTotal += fArea;

        if (writeAreas)
        {
            Info<< fArea << nl;
        }
    }

    if (writeXML)
    {
        if (writeAreas)
        {
            Info<<"</areas>" << nl;
        }

        Info<< "<area>" << areaTotal << "</area>" << nl
            << "</surfaceMeshInfo>" << nl;
    }
    else
    {
        Info<< "area      : " << areaTotal << nl;
    }

    return 0;
}

// ************************************************************************* //
