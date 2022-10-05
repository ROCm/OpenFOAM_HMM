/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "loadOrCreateMesh.H"
#include "faMesh.H"
#include "Pstream.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

Foam::boolList Foam::haveMeshFile
(
    const Time& runTime,
    const fileName& meshPath,
    const word& meshFile,
    const bool verbose
)
{
    boolList haveFileOnProc
    (
        UPstream::listGatherValues<bool>
        (
            fileHandler().isFile
            (
                fileHandler().filePath
                (
                    runTime.path()/meshPath/meshFile
                )
            )
        )
    );

    if (verbose)
    {
        Info<< "Per processor availability of \""
            << meshFile << "\" file in " << meshPath << nl
            << "    " << flatOutput(haveFileOnProc) << nl << endl;
    }

    Pstream::broadcast(haveFileOnProc);
    return haveFileOnProc;
}


void Foam::removeProcAddressing(const faMesh& mesh)
{
    IOobject io
    (
        "procAddressing",
        mesh.facesInstance(),
        faMesh::meshSubDir,
        mesh.thisDb()
    );

    for (const auto prefix : {"boundary", "edge", "face", "point"})
    {
        io.rename(prefix + word("ProcAddressing"));

        const fileName procFile(io.objectPath());
        Foam::rm(procFile);
    }
}


void Foam::removeProcAddressing(const polyMesh& mesh)
{
    IOobject io
    (
        "procAddressing",
        mesh.facesInstance(),
        polyMesh::meshSubDir,
        mesh.thisDb()
    );

    for (const auto prefix : {"boundary", "cell", "face", "point"})
    {
        io.rename(prefix + word("ProcAddressing"));

        const fileName procFile(io.objectPath());
        Foam::rm(procFile);
    }
}


void Foam::removeEmptyDir(const fileName& path)
{
    // Remove directory: silent, emptyOnly
    Foam::rmDir(path, true, true);
}


void Foam::removeEmptyDirs(const fileName& meshPath)
{
    // Delete resulting directory if empty
    fileName path(meshPath);
    path.clean();

    // Do subdirectories
    {
        const fileNameList dirs
        (
            Foam::readDir
            (
                path,
                fileName::DIRECTORY,
                false,                  // filterGz
                false                   // followLink
            )
        );
        for (const auto& dir : dirs)
        {
            removeEmptyDirs(path/dir);
        }
    }

    removeEmptyDir(path);
}


// ************************************************************************* //
