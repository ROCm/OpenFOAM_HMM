/*---------------------------------------------------------------------------*\
 =========                   |
 \\      /   F ield          | OpenFOAM: The Open Source CFD Toolbox
  \\    /    O peration      |
   \\  /     A nd            | Copyright (C) 2011-2012 OpenFOAM Foundation
    \\/      M anipulation   |
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

    As a special exception, you have permission to link this program with the
    CGAL library and distribute executables, as long as you follow the
    requirements of the GNU GPL in regard to all of the software in the
    executable aside from CGAL.

Application
    cvMesh

Description
    Conformal Voronoi automatic mesh generator

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "conformalVoronoiMesh.H"
#include "vtkSetWriter.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addBoolOption
    (
        "noFilter",
        "Do not filter the mesh"
    );
    Foam::argList::addBoolOption
    (
        "checkGeometry",
        "check all surface geometry for quality"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    runTime.functionObjects().off();

    const bool noFilter = !args.optionFound("noFilter");
    const bool checkGeometry = args.optionFound("checkGeometry");

    Info<< "Mesh filtering is " << (noFilter ? "on" : "off") << endl;

    IOdictionary cvMeshDict
    (
        IOobject
        (
            "cvMeshDict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    conformalVoronoiMesh::debug = true;

    conformalVoronoiMesh mesh(runTime, cvMeshDict);


    if (checkGeometry)
    {
        const searchableSurfaces& allGeometry = mesh.allGeometry();

        // Write some stats
        allGeometry.writeStats(List<wordList>(0), Info);
        // Check topology
        allGeometry.checkTopology(true);
        // Check geometry
        allGeometry.checkGeometry
        (
            100.0,      // max size ratio
            1e-9,       // intersection tolerance
            autoPtr<writer<scalar> >(new vtkSetWriter<scalar>()),
            0.01,       // min triangle quality
            true
        );

        return 0;
    }


    while (runTime.loop())
    {
        Info<< nl << "Time = " << runTime.timeName() << endl;

        mesh.move();

        Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    mesh.writeMesh(runTime.constant(), noFilter);

    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
