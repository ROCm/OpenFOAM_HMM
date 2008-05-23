/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "functionEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineMemberFunctionSelectionTable
    (
        functionEntry,
        insert,
        primitiveEntryIstream
    );

    defineMemberFunctionSelectionTable
    (
        functionEntry,
        insert,
        dictionaryIstream
    );
}


// * * * * * * * * * * * * Member Function Selectors * * * * * * * * * * * * //

bool Foam::functionEntry::insert
(
    const word& functionName,
    const dictionary& parentDict,
    primitiveEntry& entry,
    Istream& is
)
{
    is.fatalCheck
    (
        "functionEntry::insert"
        "(const word& functionName, const dictionary& parentDict, "
        "primitiveEntry& entry, Istream& is)"
    );

    if (!insertprimitiveEntryIstreamMemberFunctionTablePtr_)
    {
        cerr<<"functionEntry::insert"
            << "(const word&, dictionary&, primitiveEntry&, Istream&)"
            << " not yet initialized, function = "
            << functionName.c_str() << std::endl;

        // return true to keep reading anyhow
        return true;
    }

    insertprimitiveEntryIstreamMemberFunctionTable::iterator mfIter =
        insertprimitiveEntryIstreamMemberFunctionTablePtr_->find(functionName);

    if (mfIter == insertprimitiveEntryIstreamMemberFunctionTablePtr_->end())
    {
        FatalErrorIn
        (
            "functionEntry::insert"
            "(const word& functionName, const dictionary& parentDict, "
            "primitiveEntry& entry, Istream& is)"
        )   << "Unknown functionEntry " << functionName
            << endl << endl
            << "Valid functionEntries are :" << endl
            << insertprimitiveEntryIstreamMemberFunctionTablePtr_->toc()
            << exit(FatalError);
    }

    return mfIter()(parentDict, entry, is);
}


bool Foam::functionEntry::insert
(
    const word& functionName,
    dictionary& parentDict,
    Istream& is
)
{
    is.fatalCheck
    (
        "functionEntry::insert"
        "(const word& functionName, dictionary& parentDict, Istream& is)"
    );

    if (!insertdictionaryIstreamMemberFunctionTablePtr_)
    {
        cerr<<"functionEntry::insert"
            << "(const word&, dictionary&, Istream&)"
            << " not yet initialized, function = "
            << functionName.c_str() << std::endl;

        // Return true to keep reading
        return true;
    }

    insertdictionaryIstreamMemberFunctionTable::iterator mfIter =
        insertdictionaryIstreamMemberFunctionTablePtr_->find(functionName);

    if (mfIter == insertdictionaryIstreamMemberFunctionTablePtr_->end())
    {
        FatalErrorIn
        (
            "functionEntry::insert"
            "(const word& functionName, dictionary& parentDict, Istream& is)"
        )   << "Unknown functionEntry " << functionName
            << endl << endl
            << "Valid functionEntries are :" << endl
            << insertdictionaryIstreamMemberFunctionTablePtr_->toc()
            << exit(FatalError);
    }

    return mfIter()(parentDict, is);
}


// ************************************************************************* //
