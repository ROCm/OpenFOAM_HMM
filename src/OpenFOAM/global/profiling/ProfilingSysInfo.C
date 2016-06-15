/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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

#include "ProfilingSysInfo.H"
#include "demandDrivenData.H"
#include "foamVersion.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

// file-scope function
template<class T>
inline static void writeEntry
(
    Foam::Ostream& os, const Foam::word& key, const T& value
)
{
    os.writeKeyword(key) << value << Foam::token::END_STATEMENT << '\n';
}


// file-scope function
inline static void printEnv
(
    Foam::Ostream& os, const Foam::word& key, const Foam::word& envName
)
{
    const std::string value = getEnv(envName);
    if (!value.empty())
    {
        writeEntry(os, key, value);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Profiling::sysInfo::sysInfo()
{}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::Profiling::sysInfo::~sysInfo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::Ostream& Foam::Profiling::sysInfo::write
(
    Ostream& os
) const
{
    writeEntry(os, "host", hostName(false)); // short name
    writeEntry(os, "date", clock::dateTime());

    // compile-time information
    writeEntry(os, "version",     std::string(FOAMversion));
    writeEntry(os, "build",       std::string(FOAMbuild));

    printEnv(os, "arch",         "WM_ARCH");
    printEnv(os, "compilerType", "WM_COMPILER_TYPE");
    printEnv(os, "compiler",     "WM_COMPILER");
    printEnv(os, "mplib",        "WM_MPLIB");
    printEnv(os, "options",      "WM_OPTIONS");

    return os;
}


// ************************************************************************* //
