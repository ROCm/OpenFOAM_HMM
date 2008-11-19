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

Description

\*---------------------------------------------------------------------------*/

#include "SortableList.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    labelList orig(5);
    orig[0] = 7;
    orig[1] = 4;
    orig[2] = 1;
    orig[3] = 2;
    orig[4] = 9;

    labelList a(orig);
    Info << "before: " << a << endl;
    sort(a);
    Info << "after:  " << a << endl;

    SortableList<label> b(orig);
    Info << "sorted:  " << b << endl;
    Info << "indices: " << b.indices() << endl;

    Info << "shrunk:  " << b.shrink() << endl;
    Info << "indices: " << b.indices() << endl;

    // repeat by assignment
    b = orig;
    Info << "unsorted: " << b << endl;
    b.sort();

    Info << "sorted:   " << b << endl;
    Info << "indices:  " << b.indices() << endl;

    // transfer assignment
    b.transfer(orig);
    Info << "unsorted: " << b << endl;
    b.sort();

    Info << "sorted:   " << b << endl;
    Info << "indices:  " << b.indices() << endl;

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
