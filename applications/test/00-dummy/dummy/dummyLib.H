/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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

Class
    Foam::Detail::dummyLib

Description
    This class serves little useful purpose, but can be used to compile
    into a dummy library and link against for a minimal application
    that uses the wmake framework

SourceFiles
    dummyLib.C

\*---------------------------------------------------------------------------*/

#ifndef dummyLib_H
#define dummyLib_H

#include <string>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Detail
{

/*---------------------------------------------------------------------------*\
                          Class dummyLib Declaration
\*---------------------------------------------------------------------------*/

struct dummyLib
{
    //- Compile-time value of WM_ARCH
    static const std::string arch;

    //- Compile-time value of WM_COMPILER
    static const std::string compiler;

    //- Compile-time equivalent to WM_PRECISION_OPTION
    static const std::string scalar_size;

    //- Compile-time of WM_LABEL_SIZE
    static const std::string label_size;

    //- Compile-time value of WM_PRECISION_OPTION
    static const std::string precision;

    //- Compile-time value of WM_ARCH + WM_COMPILER
    static const std::string archComp;

    //- Compile-time value of WM_ARCH + WM_COMPILER, precision, label-size
    static const std::string archCompBase;

    //- DIY Compile-time value of WM_OPTIONS
    static const std::string archCompFull;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Detail
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
