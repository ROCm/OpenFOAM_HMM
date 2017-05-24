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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class UnaryMatchPredicate, class StringType>
Foam::labelList Foam::findMatchingStrings
(
    const UnaryMatchPredicate& matcher,
    const UList<StringType>& lst,
    const bool invert
)
{
    labelList indices(lst.size());

    label count = 0;
    forAll(lst, elemi)
    {
        if (matcher(lst[elemi]) ? !invert : invert)
        {
            indices[count++] = elemi;
        }
    }
    indices.setSize(count);

    return indices;
}


template<class UnaryMatchPredicate, class StringListType>
StringListType Foam::subsetMatchingStrings
(
    const UnaryMatchPredicate& matcher,
    const StringListType& lst,
    const bool invert
)
{
    // Create as a copy
    StringListType newLst(lst.size());

    // Ensure consistent addressable size (eg, DynamicList)
    newLst.setSize(lst.size());

    label count = 0;
    forAll(lst, elemi)
    {
        if (matcher(lst[elemi]) ? !invert : invert)
        {
            newLst[count++] = lst[elemi];
        }
    }
    newLst.setSize(count);

    return newLst;
}


template<class UnaryMatchPredicate, class StringListType>
void Foam::inplaceSubsetMatchingStrings
(
    const UnaryMatchPredicate& matcher,
    StringListType& lst,
    const bool invert
)
{
    label count = 0;
    forAll(lst, elemi)
    {
        if (matcher(lst[elemi]) ? !invert : invert)
        {
            if (count != elemi)
            {
                lst[count] = lst[elemi];
            }
            ++count;
        }
    }
    lst.setSize(count);
}


// ************************************************************************* //
