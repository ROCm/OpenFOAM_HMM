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
    Foam::dummyFvMesh

Description
    Functions to generate dummy finite volume meshes

SourceFiles
    dummyFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef dummyFvMesh_H
#define dummyFvMesh_H

#include "className.H"
#include "autoPtr.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;
class Time;
template<class Type> class List;

/*---------------------------------------------------------------------------*\
                         Class dummyFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dummyFvMesh
{
    //- Helper function to see if the patch type exists in the run-time
    //- selection tables
    static bool fvPatchFieldExists(const word& patchType);

    //- Set the patch information. Returns true if setting based on a real mesh
    static bool setPatchEntries
    (
        const Time& runTime,
        const word& instance,
        dictionary& patchEntries,
        label& nPatchWithFace
    );

    //- Helper function to initialise empty zones
    template<class ZoneMeshType>
    static void initialiseZone
    (
        const word& zoneTypeName,
        const fileName& instance,
        ZoneMeshType& zoneType
    );


public:

    ClassName("dummyFvMesh");

    //- Create a single cell mesh with empty boundary conditions all-round
    static autoPtr<fvMesh> singleCellMesh
    (
        const Time& runTime,
        const scalar d = 1
    );

    //- Create an equivalent 1D mesh using the same boundary types as
    //- described in the polyMesh/boundary file, complete with (empty)
    //- [point|face|cell] zones
    static autoPtr<fvMesh> equivalent1DMesh(const Time& runTime);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "dummyFvMeshTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
