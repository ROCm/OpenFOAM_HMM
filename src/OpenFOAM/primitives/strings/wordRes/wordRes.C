/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "wordRes.H"
#include "HashSet.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::wordReList Foam::wordRes::uniq(const UList<wordRe>& input)
{
    wordReList retain(input.size());
    wordHashSet uniqWord;

    label count = 0;
    forAll(input, i)
    {
        const auto& select = input[i];

        if (select.isPattern() || uniqWord.insert(select))
        {
            retain[count] = select;
            ++count;
        }
    }

    retain.setSize(count);
    return retain;
}


void Foam::wordRes::inplaceUniq(wordReList& input)
{
    wordHashSet uniqWord;

    label count = 0;
    forAll(input, i)
    {
        const auto& select = input[i];

        if (select.isPattern() || uniqWord.insert(select))
        {
            if (count != i)
            {
                input[count] = input[i];
            }
            ++count;
        }
    }

    input.setSize(count);
}


// ************************************************************************* //
