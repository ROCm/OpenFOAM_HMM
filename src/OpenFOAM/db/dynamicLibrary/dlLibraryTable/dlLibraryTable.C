/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2019 OpenCFD Ltd.
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
#include "IOstreams.H"
#include "int.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dlLibraryTable, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void* Foam::dlLibraryTable::openLibrary
(
    const fileName& libName,
    bool verbose
)
{
    if (libName.empty())
    {
        return nullptr;
    }

    std::string msg;
    void* ptr = Foam::dlOpen(fileName(libName).expand(), msg);

    DebugInFunction
        << "Opened " << libName
        << " resulting in handle " << Foam::name(ptr) << nl;

    if (!ptr)
    {
        // Even with details turned off, we want some feedback about failure
        OSstream& os = (verbose ? WarningInFunction : Serr);
        os << "Could not load " << libName << nl << msg.c_str() << endl;
    }

    return ptr;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dlLibraryTable::dlLibraryTable
(
    const UList<fileName>& libNames,
    bool verbose
)
{
    dlLibraryTable::open(libNames, verbose);
}


Foam::dlLibraryTable::dlLibraryTable
(
    const dictionary& dict,
    const word& libsEntry
)
{
    dlLibraryTable::open(dict, libsEntry);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dlLibraryTable::~dlLibraryTable()
{
    clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dlLibraryTable::empty() const
{
    for (const void* ptr : libPtrs_)
    {
        if (ptr != nullptr)
        {
            return false;
        }
    }

    return true;
}


Foam::label Foam::dlLibraryTable::size() const
{
    label nLoaded = 0;

    for (const void* ptr : libPtrs_)
    {
        if (ptr != nullptr)
        {
            ++nLoaded;
        }
    }

    return nLoaded;
}


void Foam::dlLibraryTable::clear(bool verbose)
{
    label nLoaded = 0;

    forAllReverse(libPtrs_, i)
    {
        void* ptr = libPtrs_[i];

        if (ptr == nullptr)
        {
            libNames_[i].clear();
            continue;
        }

        if (Foam::dlClose(ptr))
        {
            DebugInFunction
                << "Closed [" << i << "] " << libNames_[i]
                << " with handle " << Foam::name(ptr) << nl;

            libPtrs_[i] = nullptr;
            libNames_[i].clear();
        }
        else
        {
            ++nLoaded;  // Still loaded

            if (verbose)
            {
                WarningInFunction
                    << "Failed closing " << libNames_[i]
                    << " with handle " << Foam::name(ptr) << endl;
            }
        }
    }


    // Compact the lists
    if (nLoaded && nLoaded != libPtrs_.size())
    {
        nLoaded = 0;

        forAll(libPtrs_, i)
        {
            if (libPtrs_[i] != nullptr)
            {
                if (nLoaded != i)
                {
                    libPtrs_[nLoaded] = libPtrs_[i];
                    libNames_[nLoaded] = std::move(libNames_[i]);
                }

                ++nLoaded;
            }
        }
    }

    libPtrs_.resize(nLoaded);
    libNames_.resize(nLoaded);
}


bool Foam::dlLibraryTable::append(const fileName& libName)
{
    if (libName.empty() || libNames_.found(libName))
    {
        return false;
    }

    libPtrs_.append(nullptr);
    libNames_.append(libName);

    return true;
}


Foam::label Foam::dlLibraryTable::append(const UList<fileName>& libNames)
{
    label nAdded = 0;

    for (const fileName& libName : libNames)
    {
        if (append(libName))
        {
            ++nAdded;
        }
    }

    return nAdded;
}


bool Foam::dlLibraryTable::open(bool verbose)
{
    label nOpen = 0;
    label nCand = 0;  // Number of candidates (have libName but no pointer)

    forAll(libPtrs_, i)
    {
        const fileName& libName = libNames_[i];

        if (libPtrs_[i] == nullptr && !libName.empty())
        {
            ++nCand;
            void* ptr = openLibrary(libName, verbose);

            if (ptr)
            {
                ++nOpen;
                libPtrs_[i] = ptr;
            }
            else
            {
                libNames_[i].clear();  // Avoid trying again
            }
        }
    }

    return nOpen && nOpen == nCand;
}


void* Foam::dlLibraryTable::open
(
    const fileName& libName,
    bool verbose
)
{
    void* ptr = openLibrary(libName, verbose);

    if (ptr)
    {
        libPtrs_.append(ptr);
        libNames_.append(libName);
    }

    return ptr;
}


bool Foam::dlLibraryTable::open
(
    const UList<fileName>& libNames,
    bool verbose
)
{
    label nOpen = 0;

    for (const fileName& libName : libNames)
    {
        const label index = libNames_.find(libName);

        if (index >= 0 && libPtrs_[index] != nullptr)
        {
            // Already known and opened
            ++nOpen;
        }
        else if (dlLibraryTable::open(libName, verbose))
        {
            ++nOpen;
        }
    }

    return nOpen && nOpen == libNames.size();
}


bool Foam::dlLibraryTable::close
(
    const fileName& libName,
    bool verbose
)
{
    const label index = libNames_.rfind(libName);

    if (index < 0)
    {
        return false;
    }

    DebugInFunction
        << "Closing " << libName
        << " with handle " << Foam::name(libPtrs_[index]) << nl;

    const bool ok = Foam::dlClose(libPtrs_[index]);

    libPtrs_[index] = nullptr;
    libNames_[index].clear();

    if (!ok && verbose)
    {
        WarningInFunction
            << "Could not close " << libName << endl;
    }

    return ok;
}


void* Foam::dlLibraryTable::findLibrary(const fileName& libName)
{
    const label index = libNames_.rfind(libName);

    if (index < 0)
    {
        return nullptr;
    }

    return libPtrs_[index];
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
        if (dlLibraryTable::open(libName))  // verbose = true
        {
            ++nOpen;
        }
    }

    return nOpen && nOpen == libNames.size();
}


// ************************************************************************* //
