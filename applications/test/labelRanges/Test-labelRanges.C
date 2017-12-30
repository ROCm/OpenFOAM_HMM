/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

Description
    Test label ranges
\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "labelRanges.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::noFunctionObjects();
    argList::addArgument("start size .. startN sizeN");
    argList::addOption("verbose");
    argList::addNote
    (
        "The default is to add ranges, use 'add' and 'del' to toggle\n\n"
        "Eg, 0 10 30 10 del 20 15"
    );

    argList args(argc, argv, false, true);

    if (args.optionFound("verbose"))
    {
        labelRange::debug = 1;
    }

    {
        Info<<"test sorting" << endl;
        DynamicList<labelRange> list1(10);
        list1.append(labelRange(25, 8));
        list1.append(labelRange::identity(8));
        list1.append(labelRange(15, 5));
        list1.append(labelRange(50, -10));

        sort(list1);
        Info<<"sorted" << list1 << endl;
    }

    {
        Info<<"test intersections" << endl;
        labelRange range1(-15, 25);
        labelRange range2(7, 8);
        labelRange range3(-20, 8);
        labelRange range4(50, 8);

        Info<<range1 << " & " << range2
            << " = " << range1.subset(range2) << nl;

        Info<< range1 << " & " << range3
            << " = " << range1.subset(range3) << nl;

        Info<< range2 << " & " << range4
            << " = " << range2.subset(range4) << nl;
    }

    labelRange range;
    labelRanges ranges;

    bool removeMode = false;
    for (label argI=1; argI < args.size()-1; ++argI)
    {
        if (args[argI] == "add")
        {
            removeMode = false;
            continue;
        }
        else if (args[argI] == "del")
        {
            removeMode = true;
            continue;
        }

        {
            label start = args.argRead<label>(argI);
            label size  = args.argRead<label>(argI+1);
            ++argI;

            range.reset(start, size);
        }

        Info<< "---------------" << nl;
        if (removeMode)
        {
            Info<< "del " << range << " :";
            for (auto i : range)
            {
                Info<< " " << i;
            }
            Info<< nl;

            ranges.remove(range);
        }
        else
        {
            Info<< "add " << range  << " :";
            for (auto i : range)
            {
                Info<< " " << i;
            }
            Info<< nl;

            ranges.add(range);
        }

        Info<< "<list>" << ranges << "</list>" << nl
            << "content:";
        for (auto i : ranges)
        {
            Info<< " " << i;
        }
        Info<< nl;
    }

    return 0;
}

// ************************************************************************* //
