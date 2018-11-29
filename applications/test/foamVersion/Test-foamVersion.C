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

#include <typeinfo>
#include "foamVersion.H"
#include "Switch.H"
#include "IOstreams.H"

using namespace Foam;

int main()
{
    Info
        << "\nVersion information (function)" << nl;
    foamVersion::printBuildInfo();

    Info
        << "\nVersion information (macros)" << nl
        << "version   " << Foam::FOAMversion << nl
        << "build     " << Foam::FOAMbuild << nl
        << "buildArch " << Foam::FOAMbuildArch << nl;

    Info
        << "\nVersion information (namespace)" << nl
        << "patched?  = " << Switch(foamVersion::patched()) << nl
        << "api       " << foamVersion::api << nl
        << "patch     " << foamVersion::patch << nl
        << "version   " << foamVersion::version << nl
        << "build     " << foamVersion::build << nl
        << "buildArch " << foamVersion::buildArch << nl;

    Info
        << "\nTypes" << nl
        << "version   " << typeid(foamVersion::version).name() << nl
        << "build     " << typeid(foamVersion::build).name() << nl
        << "buildArch " << typeid(foamVersion::buildArch).name() << nl
        << "FOAMversion " << typeid(Foam::FOAMversion).name() << nl
        << "FOAMbuild   " << typeid(Foam::FOAMbuild).name() << nl;

    Info
        << "\nVerify memory addesses are identical:" << nl
        << "macro     " << long(Foam::FOAMversion) << nl
        << "namespace " << long(&(foamVersion::version[0])) << nl;

    Info
        << "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
