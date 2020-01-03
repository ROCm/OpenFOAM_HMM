/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
    Copyright (C) 2019 OpenCFD Ltd.
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
    writeMorpherCPs

Description
    Writes the NURBS3DVolume control points corresponding to dynamicMeshDict
    to csv files. Parametric coordinates are not computed.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "NURBS3DVolume.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    const dictionary NURBSdict
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                mesh.time().constant(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false
            )
        ).subDict("volumetricBSplinesMotionSolverCoeffs")
    );
    // Read box names and allocate size
    wordList controlBoxes(NURBSdict.toc());

    for (const word& boxName : controlBoxes)
    {
        if (NURBSdict.isDict(boxName))
        {
            // Creating an object writes the control points in the
            // constructor
            NURBS3DVolume::New
            (
                NURBSdict.subDict(boxName),
                mesh,
                false // do not compute parametric coordinates
            );
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
