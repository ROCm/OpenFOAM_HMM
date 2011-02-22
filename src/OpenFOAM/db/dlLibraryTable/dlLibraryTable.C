/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2011 OpenCFD Ltd.
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

#include "dlLibraryTable.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::dlLibraryTable Foam::dlLibraryTable::loadedLibraries;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dlLibraryTable::dlLibraryTable()
:
    HashTable<fileName, void*, Hash<void*> >()
{}


Foam::dlLibraryTable::readDlLibrary::readDlLibrary
(
    const dictionary& dict,
    const word& libsEntry
)
{
    open(dict, libsEntry);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dlLibraryTable::~dlLibraryTable()
{
    forAllConstIter(dlLibraryTable, *this, iter)
    {
        dlClose(iter.key());
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dlLibraryTable::open
(
    const fileName& functionLibName,
    const bool verbose
)
{
    if (functionLibName.size())
    {
        void* functionLibPtr = dlOpen(functionLibName);

        if (!functionLibPtr)
        {
            if (verbose)
            {
                WarningIn
                (
                    "dlLibraryTable::open(const fileName&)"
                )   << "could not load " << functionLibName
                    << endl;
            }

            return false;
        }
        else
        {
            if (loadedLibraries.insert(functionLibPtr, functionLibName))
            {
                return true;
            }
            else
            {
                return false;
            }
        }
    }
    else
    {
        return false;
    }
}


bool Foam::dlLibraryTable::close
(
    const fileName& functionLibName,
    const bool verbose
)
{
    void* libPtr = findLibrary(functionLibName);
    if (libPtr)
    {
        loadedLibraries.erase(libPtr);

        if (!dlClose(libPtr))
        {
            if (verbose)
            {
                WarningIn
                (
                    "dlLibraryTable::close(const fileName&)"
                )   << "could not close " << functionLibName
                    << endl;
            }

            return false;
        }

        return true;
    }
    return false;
}


void* Foam::dlLibraryTable::findLibrary(const fileName& functionLibName)
{
    forAllConstIter(dlLibraryTable, loadedLibraries, iter)
    {
        if (iter() == functionLibName)
        {
            return iter.key();
        }
    }
    return NULL;
}


bool Foam::dlLibraryTable::open
(
    const dictionary& dict,
    const word& libsEntry
)
{
    if (dict.found(libsEntry))
    {
        fileNameList libNames(dict.lookup(libsEntry));

        bool allOpened = !libNames.empty();

        forAll(libNames, i)
        {
            allOpened = dlLibraryTable::open(libNames[i]) && allOpened;
        }

        return allOpened;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
