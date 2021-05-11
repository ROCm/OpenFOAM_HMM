/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2021 OpenCFD Ltd.
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
    makeFaMesh

Description
    A mesh generator for finiteArea mesh.

Author
    Zeljko Tukovic, FAMENA
    Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#include "objectRegistry.H"
#include "Time.H"
#include "argList.H"
#include "OSspecific.H"
#include "faMesh.H"
#include "fvMesh.H"
#include "IOdictionary.H"
#include "globalIndex.H"
#include "globalMeshData.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "A mesh generator for finiteArea mesh"
    );
    argList::addOption
    (
        "empty-patch",
        "name",
        "Specify name for a default empty patch",
        false  // An advanced option, but not enough to worry about that
    );
    argList::addOption("dict", "file", "Alternative faMeshDefinition");

    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    // Reading faMeshDefinition dictionary
    #include "findMeshDefinitionDict.H"

    // Inject/overwrite name for optional 'empty' patch
    word patchName;
    if (args.readIfPresent("empty-patch", patchName))
    {
        meshDefDict.add("emptyPatch", patchName, true);
    }

    // Creation
    faMesh areaMesh(mesh, meshDefDict);

    // Writing faMesh
    Info << "Write finite area mesh ... ";
    areaMesh.write();

    #include "decomposeFaFields.H"

    Info << "\nEnd\n" << endl;

    return 0;
}

// ************************************************************************* //
