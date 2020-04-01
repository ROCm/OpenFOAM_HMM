/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

Description

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "dictionary.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #if 0
    argList::noBanner();
    argList args(argc, argv);

    if (true)
    {
        InfoErr<< "Called with " << (args.size()-1) << " args\n";
        InfoErr<< "... some error\n";
        return 2;
    }
    #endif

    FatalError.throwExceptions();

    try
    {
        WarningInFunction << "warning 1" << endl;
        IOWarningInFunction(Serr) << "warning 2" << endl;

        dictionary dict;

        IOWarningInFunction(dict) << "warning 3" << endl;

        FatalErrorInFunction
            << "This is an error from 1" << nl
            << "Explanation to follow:" << endl;

        FatalErrorInFunction
            << "Error 2"
            << exit(FatalError);
    }
    catch (const Foam::error& err)
    {
        Serr<< "Caught Foam error " << err << nl << endl;
    }

    try
    {
        FatalErrorInFunction
            << "Error# 3"
            << exit(FatalError);
    }
    catch (const Foam::error& err)
    {
        Serr<< "Caught Foam error " << err << nl << endl;
    }

    return 0;
}


// ************************************************************************* //
