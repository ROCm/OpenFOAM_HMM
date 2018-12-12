/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017-2018 OpenCFD Ltd.
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

#include "etcFiles.H"
#include "foamVersion.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * Static Functions  * * * * * * * * * * * * * * //

//
// Some of these could be exposed too (if required),
// but are fairly special purpose.
//

namespace
{

// Return the file location mode as a string.
//
// - u : location mask 0700
// - g : location mask 0070
// - o : location mask 0007
//
static inline std::string locationToString(unsigned short location)
{
    std::string mode;

    if (location & 0700) { mode += 'u'; } // User
    if (location & 0070) { mode += 'g'; } // Group
    if (location & 0007) { mode += 'o'; } // Other
    if (mode.empty()) { mode = "???"; }

    return mode;
}


// Assign 'queried' parameter to the user resource directory.
// Return true if this directory exists.
//
// Corresponds to foamEtcFile -mode=u
// Looks for
//   - ~/.OpenFOAM
static inline bool userResourceDir(Foam::fileName& queried)
{
    #ifdef FOAM_RESOURCE_USER_CONFIG_DIRNAME
    queried = Foam::home()/FOAM_RESOURCE_USER_CONFIG_DIRNAME;
    if (Foam::isDir(queried))
    {
        // If home() fails, it will have actually queried "./.OpenFOAM"
        // instead.
        // But we would have worse problems elsewhere if that were the case.
        return true;
    }
    #elif defined FULLDEBUG
        #warning FOAM_RESOURCE_USER_CONFIG_DIRNAME \
        is undefined (was this intentional?)
    #endif

    return false;
}


// Assign 'queried' parameter to the group resource directory.
// Return true if this directory exists.
// Otherwise clears the parameter and returns false.
//
// Corresponds to foamEtcFile -mode=g
// Looks for
//   - $WM_PROJECT_SITE/etc
//   - $WM_PROJECT_DIR/site/etc
static inline bool groupResourceDir(Foam::fileName& queried)
{
    #ifdef FOAM_RESOURCE_SITE_ENVNAME
    queried = Foam::getEnv(FOAM_RESOURCE_SITE_ENVNAME)/"etc";
    if (queried.size() > 3)
    {
        return Foam::isDir(queried);
    }
    #elif defined FULLDEBUG
        #warning FOAM_RESOURCE_SITE_ENVNAME \
        is undefined (was this intentional?)
    #endif

    // Fallback when WM_PROJECT_SITE is unset

    #ifdef FOAM_RESOURCE_SITE_FALLBACK_ENVNAME
    queried = Foam::getEnv(FOAM_RESOURCE_SITE_FALLBACK_ENVNAME)/"site/etc";
    if (queried.size() > 8)
    {
        return Foam::isDir(queried);
    }
    #elif defined FULLDEBUG
        #warning FOAM_RESOURCE_SITE_FALLBACK_ENVNAME \
        is undefined (was this intentional?)
    #endif

    queried.clear();
    return false;
}


// Assign 'queried' parameter to the OpenFOAM etc/ resource directory.
// Return true if it exists.
// Otherwise clears the parameter and returns false.
//
// Corresponds to foamEtcFile -mode=o
// Looks for
//   - $WM_PROJECT_DIR/etc
static inline bool projectResourceDir(Foam::fileName& queried)
{
    queried = Foam::getEnv("WM_PROJECT_DIR")/"etc";
    if (queried.size() > 3)
    {
        return Foam::isDir(queried);
    }

    queried.clear();
    return false;
}


Foam::fileNameList searchEtc
(
    const Foam::fileName& name,
    unsigned short location,
    const bool findFirst,
    bool (*accept)(const Foam::fileName&)
)
{
    // Use foamVersion::api (instead of the OPENFOAM define) to ensure this
    // stays properly synchronized with the build information
    const Foam::fileName version(std::to_string(Foam::foamVersion::api));

    Foam::fileNameList list;
    Foam::fileName queried, candidate;

    if (!(location & 0777))
    {
        // Warn about bad location (mode) ... or make it FATAL?
        std::cerr
            << "--> FOAM Error :\n    "
            "No user/group/other location specified for 'etc' file"
            " or directory\n    '"
            << name.c_str() << "'\n\n" << std::endl;
    }


    // User resource directories
    if ((location & 0700) && userResourceDir(queried))
    {
        candidate = queried/version/name;
        if (accept(candidate))
        {
            list.append(std::move(candidate));
            if (findFirst)
            {
                return list;
            }
        }

        candidate = queried/name;
        if (accept(candidate))
        {
            list.append(std::move(candidate));
            if (findFirst)
            {
                return list;
            }
        }
    }

    // Group (site) resource directories
    if ((location & 0070) && groupResourceDir(queried))
    {
        candidate = queried/version/name;
        if (accept(candidate))
        {
            list.append(std::move(candidate));
            if (findFirst)
            {
                return list;
            }
        }

        candidate = queried/name;
        if (accept(candidate))
        {
            list.append(std::move(candidate));
            if (findFirst)
            {
                return list;
            }
        }
    }

    // Other (project) resource directory
    if ((location & 0007) && projectResourceDir(queried))
    {
        candidate = queried/name;
        if (accept(candidate))
        {
            list.append(std::move(candidate));
        }
    }

    return list;
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::fileNameList Foam::etcDirs(bool test)
{
    // Use foamVersion::api (instead of the OPENFOAM define) to ensure this
    // stays properly synchronized with the build information
    const Foam::fileName version(std::to_string(Foam::foamVersion::api));

    Foam::fileNameList list(5);
    Foam::fileName queried;
    label nDirs = 0;

    // User resource directories
    if (userResourceDir(queried) || (!test && queried.size()))
    {
        list[nDirs++] = queried/version;
        list[nDirs++] = queried;
    }

    // Group (site) resource directories
    if (groupResourceDir(queried) || (!test && queried.size()))
    {
        list[nDirs++] = queried/version;
        list[nDirs++] = queried;
    }

    // Other (project) resource directory
    if (projectResourceDir(queried) || (!test && queried.size()))
    {
        list[nDirs++] = queried;
    }

    list.resize(nDirs);

    return list;
}


Foam::fileNameList Foam::findEtcDirs
(
    const fileName& name,
    const bool findFirst,
    unsigned short location
)
{
    return
        searchEtc
        (
            name,
            location,
            findFirst,
            [](const fileName& f){ return Foam::isDir(f); }
        );
}


Foam::fileNameList Foam::findEtcFiles
(
    const fileName& name,
    const bool mandatory,
    const bool findFirst,
    unsigned short location
)
{
    fileNameList list;

    if (name.size())
    {
        // A file must have a name!
        list = searchEtc
        (
            name,
            location,
            findFirst,
            [](const fileName& f){ return Foam::isFile(f); }
        );
    }

    if (mandatory && list.empty())
    {
        // Abort if file is mandatory but not found.
        // Use a direct exit, since this could occur before anything is
        // setup at all.

        std::cerr
            << "--> FOAM FATAL ERROR :\n    "
            "Could not find mandatory etc file (mode="
            << locationToString(location) << ")\n    '"
            << name.c_str() << "'\n"
            << std::endl;
        ::exit(1);
    }

    return list;
}


Foam::fileName Foam::findEtcDir
(
    const fileName& name,
    unsigned short location
)
{
    fileNameList list(findEtcDirs(name, true, location));

    fileName found;

    if (list.size())
    {
        found = std::move(list.first());
    }

    return found;
}


Foam::fileName Foam::findEtcFile
(
    const fileName& name,
    const bool mandatory,
    unsigned short location
)
{
    fileNameList list(findEtcFiles(name, mandatory, true, location));

    fileName found;

    if (list.size())
    {
        found = std::move(list.first());
    }

    return found;
}


// ************************************************************************* //
