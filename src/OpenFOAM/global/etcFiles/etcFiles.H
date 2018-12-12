/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

InNamespace
    Foam

Description
    Functions to search 'etc' directories for configuration files etc.

SourceFiles
    etcFiles.C

\*---------------------------------------------------------------------------*/

#ifndef etcFiles_H
#define etcFiles_H

#include "fileNameList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- The etc search directories in the normal search order.
//
//  \param test false to disable the default test for directory existence.
//
//  \return The list of directories
fileNameList etcDirs(bool test=true);


//- Search for directories from user/group/other directories.
//  Uses search hierarchy as per findEtcFiles().
//
//  \return The list of full paths of all the matching directories or
//      an empty list if the name cannot be found.
fileNameList findEtcDirs
(
    const fileName& name,       //!< The directory to search for
    const bool findFirst=false, //!< Stop after locating the first directory
    unsigned short location=0777 //!< User/group/other location
);


//- Search for a single directory using findEtcDirs().
//
//  \return The full path name of the first directory found in the
//      search hierarchy or an empty fileName if the name cannot be found.
fileName findEtcDir
(
    const fileName& name,       //!< The directory to search for
    unsigned short location=0777 //!< User/group/other location
);

//- Search for files from user/group/other directories.
//
//  The search hierarchy corresponds to that of the \c foamEtcFile script,
//  which allows for version-specific and version-independent files.
//
//  -# \b user settings (\c location=0700)
//    - ~/.OpenFOAM/{PROJECT_API}
//    - ~/.OpenFOAM/
//  -# \b group settings (\c location=0070)
//    - $WM_PROJECT_SITE/{PROJECT_API}/etc/
//    - $WM_PROJECT_SITE/etc/
//  -# \b other (shipped) settings (\c location=0007)
//    - $WM_PROJECT_DIR/etc/
//
//  Where {PROJECT_API} corresponds to the foamVersion::api value.
//
//  \note When \c \$WM_PROJECT_SITE is unset, uses \c $WM_PROJECT_DIR/site
//
//  \return The list of full paths of all the matching files or
//      an empty list if the name cannot be found.
fileNameList findEtcFiles
(
    const fileName& name,       //!< The file to search for
    const bool mandatory=false, //!< Abort if the file cannot be found
    const bool findFirst=false, //!< Stop after locating the first directory
    unsigned short location=0777 //!< User/group/other location
);


//- Search for a single file using findEtcFiles().
//
//  \return The full path name of the first file found in the
//      search hierarchy or an empty fileName if the name cannot be found.
fileName findEtcFile
(
    const fileName& name,       //!< The file to search for
    const bool mandatory=false, //!< Abort if the file cannot be found
    unsigned short location=0777 //!< User/group/other location
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
