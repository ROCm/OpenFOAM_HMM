/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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
    Simple test of FlatOutput

\*---------------------------------------------------------------------------*/

#include "wordList.H"
#include "ListOps.H"
#include "FlatOutput.H"
#include "IOstreams.H"
#include "macros.H"

using namespace Foam;

// For testing various pre-defined formatting
#define printFlatOutput(Content, Format) \
    STRINGIFY(Format) << ": " << flatOutput(Content, FlatOutput::Format{})


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    wordList words1
    {
        "ab", "cd", "ef", "gh",
        "ij", "kl", "mn", "op",
        "qr", "st", "uv", "wx", "yz"
    };

    {
        Info<< nl
            << "regular" << nl
            << "#----------------" << nl;

        Info<< nl << "operator<< " << words1 << nl;

        Info<< nl << "writeList: ";
        words1.writeList(Info) << nl;
    }

    Info<< nl
        << "Using c++ " << int(__cplusplus) << nl;

    {
        Info<< nl
            << "flatOutput" << nl
            << "#----------------" << nl;

        Info<< nl << "operator<< " << flatOutput(words1) << nl;

        Info<< nl << "write: ";
        flatOutput(words1).write(Info) << nl;

        Info<< nl << printFlatOutput(words1, BareComma) << nl;
        Info<< nl << printFlatOutput(words1, BareSpace) << nl;

        Info<< nl << printFlatOutput(words1, BraceComma) << nl;
        Info<< nl << printFlatOutput(words1, BraceSpace) << nl;

        Info<< nl << printFlatOutput(words1, ParenComma) << nl;
        Info<< nl << printFlatOutput(words1, ParenSpace) << nl;

        Info<< nl << printFlatOutput(words1, PointyComma) << nl;
        Info<< nl << printFlatOutput(words1, PointySpace) << nl;

        Info<< nl << printFlatOutput(words1, SquareComma) << nl;
        Info<< nl << printFlatOutput(words1, SquareSpace) << nl;
    }


    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
