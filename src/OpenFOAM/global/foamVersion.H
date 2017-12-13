/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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
    Foam

Description
    Namespace for OpenFOAM

Global
    Foam::FOAMversion

Description
    OpenFOAM version number as a static string.

Global
    Foam::FOAMbuild

Description
    OpenFOAM build information as a static string

Global
    Foam::FOAMbuildArch

Description
    OpenFOAM build architecture information (endian, label/scalar sizes)
    as a static string

SourceFiles
    global.Cver

\*---------------------------------------------------------------------------*/

#ifndef foamVersion_H
#define foamVersion_H

#include <string>

// The directory name for user-resources (located in the HOME directory)
#define WM_USER_RESOURCE_DIRNAME  ".OpenFOAM"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    extern const char* const FOAMversion;
    extern const char* const FOAMbuild;
    extern const std::string FOAMbuildArch;

    //- Additional OpenFOAM modules
    namespace Module {}
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
