/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "includeEntry.H"
#include "addToMemberFunctionSelectionTable.H"
#include "stringOps.H"
#include "IFstream.H"
#include "IOstreams.H"
#include "Time.H"
#include "fileOperation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

bool Foam::functionEntries::includeEntry::log(false);


namespace Foam
{
namespace functionEntries
{
    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        includeEntry,
        execute,
        dictionaryIstream,
        include
    );

    addNamedToMemberFunctionSelectionTable
    (
        functionEntry,
        includeEntry,
        execute,
        primitiveEntryIstream,
        include
    );
}
}

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

Foam::fileName Foam::functionEntries::includeEntry::resolveFile
(
    const fileName& dir,
    const fileName& f,
    const dictionary& dict
)
{
    fileName fName(f);

    // Substitute dictionary and environment variables.
    // Allow empty substitutions.
    stringOps::inplaceExpand(fName, dict, true, true);

    if (fName.empty() || fName.isAbsolute())
    {
        return fName;
    }

    // Relative name
    return dir/fName;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::includeEntry::execute
(
    dictionary& parentDict,
    Istream& is
)
{
    const fileName rawName(is);
    const fileName fName(resolveFile(is.name().path(), rawName, parentDict));

    autoPtr<ISstream> ifsPtr(fileHandler().NewIFstream(fName));
    ISstream& ifs = ifsPtr();

    if (ifs)
    {
        if (Foam::functionEntries::includeEntry::log)
        {
            Info<< fName << endl;
        }

        // Add watch on included file
        const dictionary& top = parentDict.topDict();
        if (isA<regIOobject>(top))
        {
            regIOobject& rio = const_cast<regIOobject&>
            (
                dynamic_cast<const regIOobject&>(top)
            );
            rio.addWatch(fName);
        }

        parentDict.read(ifs);
        return true;
    }

    FatalIOErrorInFunction
    (
        is
    )   << "Cannot open include file "
        << (ifs.name().size() ? ifs.name() : rawName)
        << " while reading dictionary " << parentDict.name()
        << exit(FatalIOError);

    return false;
}


bool Foam::functionEntries::includeEntry::execute
(
    const dictionary& parentDict,
    primitiveEntry& entry,
    Istream& is
)
{
    const fileName rawName(is);
    const fileName fName(resolveFile(is.name().path(), rawName, parentDict));

    autoPtr<ISstream> ifsPtr(fileHandler().NewIFstream(fName));
    ISstream& ifs = ifsPtr();

    if (ifs)
    {
        if (Foam::functionEntries::includeEntry::log)
        {
            Info<< fName << endl;
        }

        // Add watch on included file
        const dictionary& top = parentDict.topDict();
        if (isA<regIOobject>(top))
        {
            regIOobject& rio = const_cast<regIOobject&>
            (
                dynamic_cast<const regIOobject&>(top)
            );
            rio.addWatch(fName);
        }

        entry.read(parentDict, ifs);
        return true;
    }

    FatalIOErrorInFunction
    (
        is
    )   << "Cannot open include file "
        << (ifs.name().size() ? ifs.name() : rawName)
        << " while reading dictionary " << parentDict.name()
        << exit(FatalIOError);

    return false;
}


// ************************************************************************* //
