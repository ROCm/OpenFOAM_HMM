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
    Foam::proxyFvMesh

Description
    Functions to generate proxy finite volume meshes

SourceFiles
    proxyFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef proxyFvMesh_H
#define proxyFvMesh_H

#include "runTimeSelectionTables.H"
#include "fvMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class proxyFvMesh Declaration
\*---------------------------------------------------------------------------*/

class proxyFvMesh
:
    public fvMesh
{
protected:

    //- Helper function to initialise empty zones
    template<class ZoneMeshType>
    void initialiseZone
    (
        const word& zoneTypeName,
        const fileName& instance,
        ZoneMeshType& zoneType
    );


public:

    //- Runtime type information
    TypeName("proxyFvMesh");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            proxyFvMesh,
            time,
            (
                const Time& runTime
            ),
            (runTime)
        );


    // Selectors

        //- Return a reference to the selected proxy mesh
        static autoPtr<proxyFvMesh> New
        (
            const word& modelType,
            const Time& runTime
        );


    //- Constructor
    proxyFvMesh
    (
        const IOobject& io,
        pointField&& points,
        faceList&& faces,
        labelList&& allOwner,
        labelList&& allNeighbour
    );


    // Member Functions

        //- Helper function to see if the patch type exists in the run-time
        //- selection tables
        static bool fvPatchFieldExists(const word& patchType);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "proxyFvMeshTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
