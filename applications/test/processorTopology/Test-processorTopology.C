/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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
    Test/output processor topology

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "globalMeshData.H"
#include "OFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noFunctionObjects();
    argList::addNote
    (
        "Create graph of OpenFOAM mesh connections"
    );

    #include "setRootCase.H"

    if (!Pstream::parRun())
    {
        FatalErrorInFunction
            << "Only meaningful in parallel"
            << exit(FatalError);
    }

    #include "createTime.H"
    #include "createPolyMesh.H"

    // Adjacency table
    const labelListList& connectivity =
        mesh.globalData().topology().procNeighbours();

    if (Pstream::master())
    {
        OFstream os("processorTopology.dot");

        os << "// processorTopology" << nl << nl;
        os.beginBlock("graph");

        forAll(connectivity, proci)
        {
            label nconn = 0;

            for (const label neighProci : connectivity[proci])
            {
                if (proci < neighProci)
                {
                    if (nconn++)
                    {
                        os << "  ";
                    }
                    else
                    {
                        os << indent;
                    }
                    os << proci << " -- " << neighProci;
                }
            }

            if (nconn)
            {
                os  << nl;
            }
        }

        os.endBlock();

        Info<< "Wrote processorTopology graph: "
            << runTime.relativePath(os.name()) << nl;

        Info<< nl
            << "Use neato, circo or fdp graphviz tools" << nl;
    }

    Info<< nl << "End\n" << endl;

    return 0;
}

// ************************************************************************* //
