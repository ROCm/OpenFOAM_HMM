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
    surfaceMeshCoordinateSystemTransform

Description

    Transform (scale/rotate/translate) a surface based on
    a coordinateSystem.

Usage
    - surfaceMeshCoordinateSystemTransform inputFile outputFile [OPTION]

    @param -clean \n
    Perform some surface checking/cleanup on the input surface

    @param -scale \<scale\> \n
    Specify a scaling factor for writing the files

    @param -triSurface \n
    Use triSurface library for input/output

    @param -dict \<dictionary\> \n
    Specify an alternative dictionary for coordinateSystems.

Note
    The filename extensions are used to determine the file format type.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "Time.H"

#include "MeshedSurfaces.H"
#include "UnsortedMeshedSurfaces.H"
#include "coordinateSystems.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.append("inputFile");
    argList::validArgs.append("outputFile");
    argList::validOptions.insert("scale", "scale");
    argList::validOptions.insert("unsorted", "");
    argList::validOptions.insert("from", "sourceCoordinateSystem");
    argList::validOptions.insert("to", "targetCoordinateSystem");
    argList::validOptions.insert("dict", "dictionary");

    argList args(argc, argv);
    Time runTime(args.rootPath(), args.caseName());
    const stringList& params = args.additionalArgs();

    word dictName("coordinateSystems");
    fileName dictPath(runTime.constant());
    fileName dictLocal;

    if (args.options().found("dict"))
    {
        wordList elems(fileName(args.options()["dict"]).components());
        dictName = elems[elems.size()-1];
        dictPath = elems[0];
        dictLocal = "";

        if (elems.size() == 1)
        {
            dictPath = ".";
        }
        else if (elems.size() > 2)
        {
            dictLocal = fileName(SubList<word>(elems, elems.size()-2, 1));
        }
    }

    autoPtr<coordinateSystem> fromCsys;
    autoPtr<coordinateSystem> toCsys;

    if (args.options().found("from") || args.options().found("to"))
    {
        IOobject csDictIo
        (
            dictName,
            dictPath,
            dictLocal,
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        );

        if (!csDictIo.headerOk())
        {
            FatalErrorIn(args.executable())
                << "Cannot open coordinateSystems file\n    "
                << csDictIo.objectPath() << nl
                << exit(FatalError);
        }

        coordinateSystems csLst(csDictIo);

        if (args.options().found("from"))
        {
            word csName(args.options()["from"]);

            label csId = csLst.find(csName);
            if (csId < 0)
            {
                FatalErrorIn(args.executable())
                    << "Cannot find -from " << csName << nl
                    << "available coordinateSystems: " << csLst.toc() << nl
                    << exit(FatalError);
            }

            fromCsys.reset(new coordinateSystem(csLst[csId]));
        }

        if (args.options().found("to"))
        {
            word csName(args.options()["to"]);

            label csId = csLst.find(csName);
            if (csId < 0)
            {
                FatalErrorIn(args.executable())
                    << "Cannot find -to " << csName << nl
                    << "available coordinateSystems: " << csLst.toc() << nl
                    << exit(FatalError);
            }

            toCsys.reset(new coordinateSystem(csLst[csId]));
        }


        // maybe fix this later
        if (fromCsys.valid() && toCsys.valid())
        {
            FatalErrorIn(args.executable())
                << "Only allowed  -from  or  -to  option at the moment."
                << exit(FatalError);
        }
    }

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


    {
        meshedSurface surf(importName);

        if (args.options().found("clean"))
        {
            surf.cleanup(true);
            surf.checkOrientation(true);
        }

        if (fromCsys.valid())
        {
            tmp<pointField> tpf = fromCsys().localPosition(surf.points());
            surf.movePoints(tpf());
        }

        if (toCsys.valid())
        {
            tmp<pointField> tpf = toCsys().globalPosition(surf.points());
            surf.movePoints(tpf());
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
