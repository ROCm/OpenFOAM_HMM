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

Application
    Test-foamVersion

Description
    Print the OpenFOAM version information.

\*---------------------------------------------------------------------------*/

#include <iostream>
#include "foamVersion.H"

using namespace Foam;

int main()
{
    std::cout
        << "\nVersion information (macros)\n"
        << "version   " << Foam::FOAMversion << '\n'
        << "build     " << Foam::FOAMbuild << '\n'
        << "buildArch " << Foam::FOAMbuildArch << '\n';

    std::cout
        << "\nVersion information (namespace)\n"
        << "version   " << foamVersion::version << '\n'
        << "build     " << foamVersion::build << '\n'
        << "buildArch " << foamVersion::buildArch << '\n';

    std::cout
        << "\nVerify memory addesses are identical:\n"
        << "macro     " << long(&(Foam::FOAMversion)) << '\n'
        << "namespace " << long(&(foamVersion::version)) << '\n';

    std::cout
        << "\nEnd\n";

    return 0;
}

// ************************************************************************* //
