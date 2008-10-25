/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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

Application
    surfaceMeshConvert

Description
    Converts from one surface mesh format to another

Usage
    - surfaceMeshConvert inputFile outputFile [OPTION]

    @param -clean \n
    Perform some surface checking/cleanup on the input surface

    @param -scale \<scale\> \n
    Specify a scaling factor for writing the files

    @param -triSurface \n
    Use triSurface library for input/output

Note
    The filename extensions are used to determine the file format type.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "Time.H"
#include "polyMesh.H"
#include "meshedSurface.H"
#include "triSurface.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.append("inputFile");
    argList::validArgs.append("outputFile");
    argList::validOptions.insert("clean", "");
    argList::validOptions.insert("scale", "scale");
    argList::validOptions.insert("triSurface", "");
#   include "setRootCase.H"
    const stringList& params = args.additionalArgs();

    scalar scaleFactor = 0;
    if (args.options().found("scale"))
    {
        IStringStream(args.options()["scale"])() >> scaleFactor;
    }

    fileName importName(params[0]);
    fileName exportName(params[1]);

    if (importName == exportName)
    {
        FatalErrorIn(args.executable())
            << "Output file " << exportName << " would overwrite input file."
            << exit(FatalError);
    }

    if
    (
        !meshedSurface::canRead(importName.ext(), true)
     || !meshedSurface::canWrite(exportName.ext(), true)
    )
    {
        return 1;
    }

    if (args.options().found("triSurface"))
    {
        triSurface surf(importName);

        if (args.options().found("clean"))
        {
            surf.cleanup(true);
            surf.checkOrientation(true);
        }

        Info << "writing " << exportName;
        if (scaleFactor <= 0)
        {
            Info<< " without scaling" << endl;
        }
        else
        {
            Info<< " with scaling " << scaleFactor << endl;
            surf.scalePoints(scaleFactor);
        }
        surf.write(exportName);
    }
    else
    {
        meshedSurface surf(importName);

        if (args.options().found("clean"))
        {
            surf.cleanup(true);
            surf.checkOrientation(true);
        }

        Info<< "writing " << exportName;
        if (scaleFactor <= 0)
        {
            Info<< " without scaling" << endl;
        }
        else
        {
            Info<< " with scaling " << scaleFactor << endl;
            surf.scalePoints(scaleFactor);
        }
        surf.write(exportName);
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}

// ************************************************************************* //
