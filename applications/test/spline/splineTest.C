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
#include "argList.H"

#include "vector.H"
#include "IFstream.H"
#include "BSpline.H"
#include "CatmullRomSpline.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.insert("file .. fileN");
    argList::addBoolOption("B", "B-Spline");
    argList::addBoolOption("cmr", "catmull-rom spline (default)");
    argList::addOption
    (
        "n",
        "INT",
        "number of segments for evaluation - default 20"
    );

    argList args(argc, argv, false, true);

    if (args.additionalArgs().empty())
    {
        args.printUsage();
    }

    bool useBSpline    = args.optionFound("B");
    bool useCatmullRom = args.optionFound("cmr");
    label nSeg = args.optionLookupOrDefault<label>("n", 20);

    if (!useBSpline && !useCatmullRom)
    {
        Info<<"defaulting to Catmull-Rom spline" << endl;
        useCatmullRom = true;
    }

    forAll(args.additionalArgs(), argI)
    {
        const string& srcFile = args.additionalArgs()[argI];
        Info<< nl << "reading " << srcFile << nl;
        IFstream ifs(srcFile);

        List<pointField> pointFields(ifs);

        forAll(pointFields, splineI)
        {
            Info<<"\noriginal points: " << pointFields[splineI] << nl;

            if (useBSpline)
            {
                BSpline spl(pointFields[splineI], vector::zero, vector::zero);

                Info<< nl
                    << "B-Spline interpolation:" << nl
                    << "----------------------" << endl;

                for (label segI = 0; segI <= nSeg; ++segI)
                {
                    scalar lambda = scalar(segI)/scalar(nSeg);
                    Info<< spl.position(lambda) << "    // " << lambda << endl;
                }
            }

            if (useCatmullRom)
            {
                CatmullRomSpline spl
                (
                    pointFields[splineI]
                );

                Info<< nl
                    <<"Catmull-Rom interpolation:" << nl
                    << "-------------------------" << endl;

                for (label segI = 0; segI <= nSeg; ++segI)
                {
                    scalar lambda = scalar(segI)/scalar(nSeg);
                    Info<< spl.position(lambda) << "    // " << lambda << endl;
                }
            }
        }
    }

    return 0;
}


// ************************************************************************* //
