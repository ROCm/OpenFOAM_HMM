/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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
    Test allocation patterns when reading into an existing list.

\*---------------------------------------------------------------------------*/

#include "DynamicList.H"
#include "DynamicField.H"
#include "IOstreams.H"
#include "ITstream.H"
#include "OTstream.H"
#include "StringStream.H"
#include "FlatOutput.H"
#include "ListOps.H"
#include "labelRange.H"
#include "labelIndList.H"

using namespace Foam;

template<class T, int SizeMin>
void printInfo
(
    const word& tag,
    const DynamicList<T, SizeMin>& list,
    const bool showSize = true
)
{
    Info<< '<' << tag;
    if (showSize)
    {
        Info<< " size=\"" << list.size()
            << "\" capacity=\"" << list.capacity() << "\"";
        if (list.cdata())
        {
            Info<< " ptr=\"" << name(list.cdata()) << "\"";
        }
        else
        {
            Info<< " ptr=\"nullptr\"";
        }
    }
    Info<< '>' << nl << flatOutput(list) << nl
        << "</" << tag << ">\n" << endl;
}


template<class T, int SizeMin>
void printInfo
(
    const word& tag,
    const DynamicField<T, SizeMin>& list,
    const bool showSize = true
)
{
    Info<< '<' << tag;
    if (showSize)
    {
        Info<< " size=\"" << list.size()
            << "\" capacity=\"" << list.capacity() << "\"";
        if (list.cdata())
        {
            Info<< " ptr=\"" << name(list.cdata()) << "\"";
        }
        else
        {
            Info<< " ptr=\"nullptr\"";
        }
    }
    Info<< '>' << nl << flatOutput(list) << nl
        << "</" << tag << ">\n" << endl;
}


template<class T, int SizeMin>
void readList
(
    DynamicList<T, SizeMin>& output,
    const UList<T>& input
)
{
    OTstream os;
    os << input;
    ITstream is("input", os.tokens());

    is >> output;
}

template<class T, int SizeMin>
void readList
(
    DynamicField<T, SizeMin>& output,
    const UList<T>& input
)
{
    OTstream os;
    os << input;
    ITstream is("input", os.tokens());

    is >> output;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    //
    {
        DynamicList<label, 64> list1;

        list1.resize(4);
        ListOps::identity(list1);

        list1.resize(3);
        printInfo("", list1);

        // list1.clear();
        // printInfo("", list1);

        list1.setCapacity(3);
        printInfo("", list1);
    }

    Info<< "\nEnd\n";

    return 0;
}


// ************************************************************************* //
