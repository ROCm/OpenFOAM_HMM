/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2020 OpenCFD Ltd.
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
    Minimal compilation test with wmake, without OpenFOAM libraries.

    The application and library can also serve as a minimal test case for
    wmake, or to provide a minimal library/executable target for testing.

\*---------------------------------------------------------------------------*/

#include "dummyLib.H"
#include <cstring>
#include <iostream>

constexpr char nl = '\n';
constexpr const char* const bold = "\\fB";  // nroff
constexpr const char* const norm = "\\fR";  // nroff

constexpr const char* const website = "www.openfoam.com";

using std::cout;
using wmake = Foam::Detail::dummyLib;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

static void printMan(const char* exeName)
{
    cout
        << ".TH \"" << exeName << "\" 1 "
        << "\"OpenFOAM-v" << OPENFOAM
        << "\" \"" << website << "\" \"OpenFOAM Commands Manual\""
        << nl;

    cout
        << ".SH NAME" << nl
        << exeName
        << " \\- part of " << bold << "OpenFOAM" << norm
        << " (The Open Source CFD Toolbox)." << nl
        << ".SH SYNOPSIS" << nl
        << bold << exeName << norm << " [OPTIONS]" << nl;

    cout
        << ".SH DESCRIPTION" << nl
        << ".nf" << nl
        << "Minimal compilation test with wmake, without OpenFOAM libraries."
        << nl
        << ".fi" << nl;

    cout
        << ".SH OPTIONS" << nl
        << ".TP" << nl
        << "-help-man" << nl
        << "Display manpage" << nl;

    cout
        << ".SH INFORMATION" << nl
        << ".nf" << nl
        << "label    = " << wmake::label_size << nl
        << "scalar   = " << wmake::scalar_size;

    if
    (
        wmake::solveScalar_size
     && wmake::solveScalar_size != wmake::scalar_size
    )
    {
        cout
            << " [solve=" << wmake::solveScalar_size << "]";
    }
    cout
        << " (" << wmake::precision << ')' << nl
        << "arch     = " << wmake::arch << nl
        << "compiler = " << wmake::compiler << nl;

    cout
        << nl
        << "archComp     = " << wmake::archComp << nl
        << "archCompBase = " << wmake::archCompBase << nl
        << "archCompFull = " << wmake::archCompFull << nl;
    cout
        << ".fi" << nl;

    cout
        << ".SH \"SEE ALSO\"" << nl
        << "Online documentation https://" << website << "/documentation/"
        << nl;
}


int main(int argc, char *argv[])
{
    // Process -help-man
    if (argc > 1 && strcmp(argv[1], "-help-man") == 0)
    {
        printMan("Test-dummyLib");
        return 0;
    }

    cout
        << nl
        << "OPENFOAM  = " << OPENFOAM << nl
        << "label     = " << wmake::label_size << nl
        << "scalar    = " << wmake::scalar_size
        << " (" << wmake::precision << ')' << nl;

    if
    (
        wmake::solveScalar_size
     && wmake::solveScalar_size != wmake::scalar_size
    )
    {
        cout
            << "solve     = " << wmake::solveScalar_size << nl;
    }

    cout
        << "arch      = " << wmake::arch << nl
        << "compiler  = " << wmake::compiler << nl;

    cout
        << nl
        << "archComp     = " << wmake::archComp << nl
        << "archCompBase = " << wmake::archCompBase << nl
        << "archCompFull = " << wmake::archCompFull << nl;

    cout<< nl;

    return 0;
}


// ************************************************************************* //
