/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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

#include "IOobjectList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Templated implementation for classes()
template<class MatchPredicate>
Foam::HashTable<Foam::wordHashSet> Foam::IOobjectList::classesImpl
(
    const IOobjectList& list,
    const MatchPredicate& matchName
)
{
    HashTable<wordHashSet> summary(2*list.size());

    // Summary (key,val) = (class-name, object-names)
    forAllConstIters(list, iter)
    {
        const word& key = iter.key();
        const IOobject* io = iter.object();

        if (matchName(key))
        {
            // Create entry (if needed) and insert
            summary(io->headerClassName()).insert(key);
        }
    }

    return summary;
}


// Templated implementation for names(), sortedNames()
template<class MatchPredicate>
Foam::wordList Foam::IOobjectList::namesImpl
(
    const IOobjectList& list,
    const word& matchClass,
    const MatchPredicate& matchName,
    const bool doSort
)
{
    wordList objNames(list.size());

    label count = 0;
    forAllConstIters(list, iter)
    {
        const word& key = iter.key();
        const IOobject* io = iter.object();

        if (matchClass(io->headerClassName()) && matchName(key))
        {
            objNames[count] = key;
            ++count;
        }
    }

    objNames.setSize(count);

    if (doSort)
    {
        Foam::sort(objNames);
    }

    return objNames;
}


// Templated implementation for lookup()
template<class MatchPredicate>
Foam::IOobjectList Foam::IOobjectList::lookupImpl
(
    const IOobjectList& list,
    const MatchPredicate& matchName
)
{
    IOobjectList results(list.size());

    forAllConstIters(list, iter)
    {
        const word& key = iter.key();
        const IOobject* io = iter.object();

        if (matchName(key))
        {
            if (IOobject::debug)
            {
                InfoInFunction << "Found " << key << endl;
            }

            results.set(key, new IOobject(*io));
        }
    }

    return results;
}


// Templated implementation for lookupClass()
template<class MatchPredicate>
Foam::IOobjectList Foam::IOobjectList::lookupClassImpl
(
    const IOobjectList& list,
    const word& matchClass,
    const MatchPredicate& matchName
)
{
    IOobjectList results(list.size());

    forAllConstIters(list, iter)
    {
        const word& key = iter.key();
        const IOobject* io = iter.object();

        if (matchClass(io->headerClassName()) && matchName(key))
        {
            if (IOobject::debug)
            {
                InfoInFunction << "Found " << key << endl;
            }

            results.set(key, new IOobject(*io));
        }
    }

    return results;
}


// ************************************************************************* //
