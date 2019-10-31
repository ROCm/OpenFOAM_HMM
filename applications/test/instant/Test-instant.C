/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd.
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
    Test instant, fileNameInstant

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "instant.H"
#include "fileNameInstant.H"
#include "DynamicList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    DynamicList<instant> times;

    times.append(instant{});
    times.append({12, "abc"});
    times.append(instant{3.14159});
    times.append({300.456, "def"});
    times.append({454.456, "xyz"});
    times.append({10, "ten"});

    {
        word timeName("twenty");

        Info<<"move append: " << timeName << nl;
        times.append({20, std::move(timeName)});
        Info<<"after append: " << timeName << nl;
    }

    Info<< nl << "times:" << times << nl;
    sort(times);
    Info<< "Sorted:" << times << nl;


    DynamicList<fileNameInstant> files;
    files.append(fileNameInstant{});
    files.append({12, "twelve"});
    files.append({3.14, "/path/almost-pi"});
    files.append({300, "/dev/value"});
    files.append({454, "/tmp/xyz"});
    files.append({10, "ten"});

    Info<< nl << "files:" << files << nl;
    sort(files);
    Info<< "Sorted:" << files << nl;


    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
