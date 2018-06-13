/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2018 OpenCFD Ltd.
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
#include "int.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dlLibraryTable, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dlLibraryTable::dlLibraryTable()
{}


Foam::dlLibraryTable::dlLibraryTable
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
    forAllReverse(libPtrs_, i)
    {
        if (libPtrs_[i])
        {
            if (debug)
            {
                InfoInFunction
                    << "Closing " << libNames_[i]
                    << " with handle " << uintptr_t(libPtrs_[i]) << endl;
            }
            if (!dlClose(libPtrs_[i]))
            {
                WarningInFunction<< "Failed closing " << libNames_[i]
                    << " with handle " << uintptr_t(libPtrs_[i]) << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dlLibraryTable::open
(
    const fileName& libName,
    const bool verbose
)
{
    if (libName.empty())
    {
        return false;
    }

    void* ptr = dlOpen(fileName(libName).expand(), verbose);

    if (debug)
    {
        InfoInFunction
            << "Opened " << libName
            << " resulting in handle " << uintptr_t(ptr) << endl;
    }

    if (ptr)
    {
        libPtrs_.append(ptr);
        libNames_.append(libName);
        return true;
    }

    if (verbose)
    {
        WarningInFunction
            << "could not load " << libName
            << endl;
    }

    return false;
}


bool Foam::dlLibraryTable::close
(
    const fileName& libName,
    const bool verbose
)
{
    label index = -1;
    forAllReverse(libNames_, i)
    {
        if (libName == libNames_[i])
        {
            index = i;
            break;
        }
    }

    if (index == -1)
    {
        return false;
    }

    if (debug)
    {
        InfoInFunction
            << "Closing " << libName
            << " with handle " << uintptr_t(libPtrs_[index]) << endl;
    }

    const bool ok = dlClose(libPtrs_[index]);

    libPtrs_[index] = nullptr;
    libNames_[index].clear();

    if (!ok && verbose)
    {
        WarningInFunction
            << "could not close " << libName
            << endl;
    }

    return ok;
}


void* Foam::dlLibraryTable::findLibrary(const fileName& libName)
{
    label index = -1;
    forAllReverse(libNames_, i)
    {
        if (libName == libNames_[i])
        {
            index = i;
            break;
        }
    }

    if (index != -1)
    {
        return libPtrs_[index];
    }

    return nullptr;
}


bool Foam::dlLibraryTable::open
(
    const dictionary& dict,
    const word& libsEntry
)
{
    fileNameList libNames;
    dict.readIfPresent(libsEntry, libNames);

    label nOpen = 0;

    for (const fileName& libName : libNames)
    {
        if (dlLibraryTable::open(libName))
        {
            ++nOpen;
        }
    }

    return nOpen && nOpen == libNames.size();
}


// ************************************************************************* //
