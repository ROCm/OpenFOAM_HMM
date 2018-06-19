/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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
    Test-coordinateSystem

Description
    Expand coordinate system definitions

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "coordinateSystem.H"
#include "Fstream.H"
#include "IOstreams.H"

using namespace Foam;

void doTest(const dictionary& dict)
{
    Info<< dict.dictName() << dict << nl;

    // Could fail?
    const bool throwingIOError = FatalIOError.throwExceptions();
    const bool throwingError = FatalError.throwExceptions();
    try
    {
        coordinateSystem cs1(dict.dictName(), dict);

        coordinateSystem cs2;

        // Move assign
        cs2 = std::move(cs1);

        // Info<<cs2 << nl;
        cs2.writeDict(Info, true);
        Info<< nl;
    }
    catch (Foam::IOerror& err)
    {
        Info<< "Caught FatalIOError " << err << nl << endl;
    }
    catch (Foam::error& err)
    {
        Info<< "Caught FatalError " << err << nl << endl;
    }
    FatalError.throwExceptions(throwingError);
    FatalIOError.throwExceptions(throwingIOError);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addArgument("dict .. dictN");
    argList args(argc, argv, false, true);

    if (args.size() <= 1)
    {
        Info<<"no coordinateSystem dictionaries to expand" << nl;
    }
    else
    {
        for (label argi=1; argi < args.size(); ++argi)
        {
            const string& dictFile = args[argi];
            IFstream is(dictFile);

            dictionary inputDict(is);

            forAllConstIters(inputDict, iter)
            {
                if (iter().isDict())
                {
                    doTest(iter().dict());
                }
            }
        }
    }

    return 0;
}


// ************************************************************************* //
