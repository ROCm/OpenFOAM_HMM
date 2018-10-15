/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2018 OpenCFD Ltd.
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

Description

\*---------------------------------------------------------------------------*/

#include "OSspecific.H"
#include "clock.H"
#include "clockTime.H"
#include "cpuTime.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    Info<<"cwd() " << cwd() << nl;
    Info<<"cwd(-P) " << cwd(false) << nl;
    Info<<"cwd(-L) " << cwd(true) << nl;

    Info<<"rmDir" << nl;
    rmDir("hmm");

    {
        Foam::clock sysClock();

        Info<< "clock: date "  << clock::date() << nl
            << "clock: time "  << clock::clockTime() << nl
            << "clock: iso  "  << clock::dateTime() << nl;
    }

    Info<< "since epoch = " << clockValue::now().str() << nl;

    {
        clockValue a;

        Info<< "clockValue() " << a << nl;
        a.update();
        Info<< "updated " << a << nl;

        Info<< "sleep 4..." << endl;
        sleep(4);

        a.update();
        Info<< " = " << a.seconds() << nl;

        Info<< "sleep 2..." << endl;
        sleep(2);

        Info<< "elapsed = " << a.elapsed() << nl;
        Info<< "elapsed = " << a.elapsed().seconds() << nl;
        Info<< "elapsed = " << a.elapsed().str() << nl;

        clockValue b = clockValue::now();

        Info<< "(" << b << " - " << a << ") = " << (b - a) << nl;
        Info<< "(" << b << " + " << a << ") = " << (b + a) << nl;
    }

    {
        clockTime clk;

        Info<< "starting clockTime" << nl;

        Info<< "sleep 4..." << endl;
        sleep(4);

        Info<< "increment = " << clk.timeIncrement() << nl;
        Info<< "elapsed   = " << clk.elapsedTime() << nl;

        Info<< "sleep 4..." << endl;
        sleep(4);
        Info<< "increment = " << clk.timeIncrement() << nl;
        Info<< "elapsed   = " << clk.elapsedTime() << nl;

        Info<< "sleep 2..." << endl;
        sleep(2);
        Info<< "elapsed   = " << clk.elapsedTime() << nl;
        Info<< "increment = " << clk.timeIncrement() << nl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
