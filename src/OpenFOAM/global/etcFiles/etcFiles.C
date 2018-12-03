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
// These could be exposed too (if required), but are fairly special purpose.
//

namespace
{

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

    return false;
}


// Assign 'queried' parameter to the OpenFOAM etc/ resource directory.
// Return true if it exists.
//
// Corresponds to foamEtcFile -mode=o
// Looks for
//   - $WM_PROJECT_DIR/etc
static inline bool projectResourceDir(Foam::fileName& queried)
{
    queried = Foam::getEnv("WM_PROJECT_DIR")/"etc";
    return (queried.size() > 3 && Foam::isDir(queried));
}


Foam::fileNameList searchEtc
(
    const Foam::fileName& name,
    const bool findFirst,
    bool (*accept)(const Foam::fileName&)
)
{
    Foam::fileName version(Foam::getEnv("WM_PROJECT_VERSION"));

    // Fallback when WM_PROJECT_VERSION is unset
    if (version.empty())
    {
        #if OPENFOAM
        version.assign(std::to_string(OPENFOAM));
        #else
        version.assign(foamVersion::version);
        #endif
    }


    Foam::fileNameList list;
    Foam::fileName dir, candidate;


    // User resource directories
    if (userResourceDir(dir))
    {
        candidate = dir/version/name;
        if (accept(candidate))
        {
            list.append(std::move(candidate));
            if (findFirst)
            {
                return list;
            }
        }

        candidate = dir/name;
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
    if (groupResourceDir(dir))
    {
        candidate = dir/version/name;
        if (accept(candidate))
        {
            list.append(std::move(candidate));
            if (findFirst)
            {
                return list;
            }
        }

        candidate = dir/name;
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
    if (projectResourceDir(dir))
    {
        candidate = dir/name;
        if (accept(candidate))
        {
            list.append(std::move(candidate));
        }
    }

    return list;
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::fileNameList Foam::findEtcDirs
(
    const fileName& name,
    const bool findFirst
)
{
    return
        searchEtc
        (
            name,
            findFirst,
            [](const fileName& f){ return Foam::isDir(f); }
        );
}


Foam::fileNameList Foam::findEtcFiles
(
    const fileName& name,
    const bool mandatory,
    const bool findFirst
)
{
    fileNameList list;

    if (name.size())
    {
        // A file must have a name!
        list = searchEtc
        (
            name,
            findFirst,
            [](const fileName& f){ return Foam::isFile(f); }
        );
    }

    if (mandatory && list.empty())
    {
        // Abort if file is mandatory but not found
        std::cerr
            << "--> FOAM FATAL ERROR in Foam::findEtcFiles()"
               " :  could not find mandatory file\n    '"
            << name.c_str() << "'\n\n" << std::endl;
        ::exit(1);
    }

    return list;
}


Foam::fileName Foam::findEtcFile(const fileName& name, const bool mandatory)
{
    fileNameList list(findEtcFiles(name, mandatory, true));

    if (list.size())
    {
        return list.first();
    }

    return fileName();
}


// ************************************************************************* //
