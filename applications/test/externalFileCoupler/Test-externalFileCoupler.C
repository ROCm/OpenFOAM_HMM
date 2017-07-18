/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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
    Test-externalFileCoupler

Description
    Test of master/slave communication etc.
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "externalFileCoupler.H"

using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::addOption("max", "N", "max number of calls (default: 1000)");
    argList::addBoolOption("slave", "run as slave");

    #include "setRootCase.H"

    const label maxCount = args.optionLookupOrDefault<label>("max", 1000);

    externalFileCoupler coupler;

    if (args.optionFound("slave"))
    {
        const word role = "slave";
        Info<< "Running as " << role << " max=" << maxCount << endl;

        for (label count = 0; count < maxCount; ++count)
        {
            // Wait for master, but stop if status=done was seen

            Info<< role << ": waiting for master" << endl;
            if (!coupler.waitForMaster())
            {
                Info<< role << ": stopping. status=done was detected" << endl;
                break;
            }

            Info<< role << ": switch to master" << endl;
            coupler.useMaster();
        }
    }
    else
    {
        const word role = "master";
        Info<< "Running as " << role << " max=" << maxCount << endl;

        for (label count = 0; count < maxCount; ++count)
        {
            // Wait for slave, but stop if status=done was seen

            Info<< role << ": waiting for slave" << endl;
            if (!coupler.waitForSlave())
            {
                Info<< role << ": stopping. status=done was detected" << endl;
                break;
            }

            Info<< role << ": switch to slave" << endl;
            coupler.useSlave();
        }
    }

    return 0;
}


// ************************************************************************* //
