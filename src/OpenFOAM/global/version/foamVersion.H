/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Namespace
    Foam::foamVersion

Description
    Namespace for OpenFOAM version information

Note
    Compile-time version information is conveyed by the \b OPENFOAM define
    provided in the wmake rules "General/general".
    This value is mirrored by the \b OPENFOAM_COM and \b OPENFOAM_PLUS defines
    provided via the stdFoam.H file.
    The OPENFOAM_PLUS define was last used for 1712.

SourceFiles
    foamVersion.C
    global.Cver

\*---------------------------------------------------------------------------*/

#ifndef foamVersion_H
#define foamVersion_H

#include "stdFoam.H"
#include <string>

// The directory name for user-resources (located in the HOME directory)
#define WM_USER_RESOURCE_DIRNAME  ".OpenFOAM"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    //- Version information
    namespace foamVersion
    {
        //- OpenFOAM version (number or name) as a static c-string.
        extern const char* const version;

        //- OpenFOAM build information as a static c-string
        extern const char* const build;

        //- OpenFOAM build architecture information
        //- (machine endian, label/scalar sizes) as a static std::string
        extern const std::string buildArch;

        //- Print information about version, build, arch to Info
        //
        //  \verbatim
        //  Using: OpenFOAM-<VERSION> (see www.OpenFOAM.com)
        //  Build: <BUILD> (OPENFOAM=<API_LEVEL>)
        //  Arch:  "<ARCH_INFO>"
        //  \endverbatim
        //
        //  \param full includes Arch information
        void printBuildInfo(const bool full=true);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Compatibility names (1806 and earlier)

#define FOAMversion     foamVersion::version
#define FOAMbuild       foamVersion::build
#define FOAMbuildArch   foamVersion::buildArch

#endif

// ************************************************************************* //
