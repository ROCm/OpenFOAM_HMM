/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2012-2013 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "wordIOList.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineCompoundTypeName(List<word>, wordList);
    addCompoundToRunTimeSelectionTable(List<word>, wordList);

    defineTemplateTypeNameAndDebugWithName(wordIOList, "wordList", 0);
    defineTemplateTypeNameAndDebugWithName(wordListIOList, "wordListList", 0);
}


void Foam::printTable
(
    const List<wordList>& wll,
    List<string::size_type>& columnWidth,
    Ostream& os
)
{
    if (wll.empty()) return;

    // Find the maximum word length for each column
    columnWidth.setSize(wll[0].size(), string::size_type(0));
    forAll(columnWidth, coli)
    {
        forAll(wll, rowi)
        {
            columnWidth[coli] =
                std::max
                (
                    columnWidth[coli],
                    string::size_type(wll[rowi][coli].size())
                );
        }
    }

    // Print the rows adding spacing for the columns
    forAll(wll, rowi)
    {
        forAll(wll[rowi], coli)
        {
            os  << wll[rowi][coli];
            for
            (
                string::size_type space=0;
                space < columnWidth[coli] - wll[rowi][coli].size() + 2;
                ++space
            )
            {
                os  << ' ';
            }
        }
        os  << nl;

        if (!rowi) os  << nl;
    }
}


void Foam::printTable(const List<wordList>& wll, Ostream& os)
{
    List<string::size_type> columnWidth;
    printTable(wll, columnWidth, os);
}


// ************************************************************************* //
