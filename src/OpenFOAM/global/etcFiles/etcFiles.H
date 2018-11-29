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

//- Search for directories from user/group/other directories.
//  Uses search hierarchy as per findEtcFiles().
//
//  \return The list of full paths of all the matching directories or
//      an empty list if the name cannot be found.
fileNameList findEtcDirs
(
    const fileName& name=fileName::null, //!< the file to search for
    const bool findFirst=false  //!< stop when the first file has been found
);


//- Search for files from user/group/other directories.
//
//  The search hierarchy corresponds to that of the foamEtcFile script,
//  which allows for version-specific and version-independent files:
//  -# \b user settings
//    - ~/.OpenFOAM/{PROJECT_API}
//    - ~/.OpenFOAM/
//  -# \b group settings
//    - $WM_PROJECT_SITE/{PROJECT_API}
//    - $WM_PROJECT_SITE/
//  -# \b other (shipped) settings
//    - $WM_PROJECT_DIR/etc/
//
//  Where {PROJECT_API} is the value of the OPENFOAM define.
//  - \b \$WM_PROJECT_SITE : If unset, use $WM_PROJECT_DIR/site
//
//  \return The list of full paths of all the matching files or
//      an empty list if the name cannot be found.
fileNameList findEtcFiles
(
    const fileName& name,       //!< the file to search for
    const bool mandatory=false, //!< abort if the file cannot be found
    const bool findFirst=false  //!< stop when the first file has been found
);


//- Search for a single file using findEtcFiles().
//
//  \return The full path name of the first file found in the
//      search hierarchy or an empty fileName if the name cannot be found.
fileName findEtcFile
(
    const fileName& name,       //!< the file to search for
    const bool mandatory=false  //!< abort if the file cannot be found
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
