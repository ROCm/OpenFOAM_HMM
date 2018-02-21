/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2018 OpenCFD Ltd.
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

Foam::wordRes Foam::wordRes::uniq(const UList<wordRe>& input)
{
    wordRes output(input.size());
    wordHashSet uniqWord;

    label count = 0;
    for (const wordRe& select : input)
    {
        if (select.isPattern() || uniqWord.insert(select))
        {
            output[count] = select;
            ++count;
        }
    }

    output.resize(count);
    return output;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::wordRes::uniq()
{
    wordHashSet uniqWord;

    label i = 0, count = 0;
    for (wordRe& select : *this)
    {
        if (select.isPattern() || uniqWord.insert(select))
        {
            if (count != i)
            {
                (*this)[count] = std::move(select);
            }
            ++count;
        }
        ++i;
    }

    resize(count);
}


// ************************************************************************* //
