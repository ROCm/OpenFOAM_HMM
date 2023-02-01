/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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
    Test loading of different gravity items

\*---------------------------------------------------------------------------*/

#include "MeshObject.H"
#include "gravityMeshObject.H"
#include "IOobjectList.H"
#include "IOstreams.H"
#include "argList.H"
#include "Time.H"

using namespace Foam;

void printInfo(const meshObjects::gravity& g)
{
    Pout<< "name:" << g.uniformDimensionedVectorField::name()
        << " type:" << g.type()
        << " value:" << g.value() << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    IOobjectList objects(runTime, runTime.constant());

    Info<< "Found: " << objects << nl << endl;

    for (const IOobject& io : objects.sorted<uniformDimensionedVectorField>())
    {
        if (io.name() == meshObjects::gravity::typeName)
        {
            const auto& g = meshObjects::gravity::New(runTime);
            printInfo(g);
        }
        else
        {
            const auto& g = meshObjects::gravity::New(io.name(), runTime);
            printInfo(g);
        }

        Pout<< "registered:" << flatOutput(runTime.sortedToc()) << nl << endl;
    }

    meshObjects::gravity::Delete("g", runTime);
    meshObjects::gravity::Delete("something-not-in-registry", runTime);

    Info<< "after Delete" << nl;
    Pout<< "registered:" << flatOutput(runTime.sortedToc()) << endl;

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
