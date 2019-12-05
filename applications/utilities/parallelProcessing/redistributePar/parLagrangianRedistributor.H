/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::parLagrangianRedistributor

Description
    Lagrangian field redistributor.

    Runs in parallel. Redistributes from fromMesh to toMesh.

SourceFiles
    parLagrangianRedistributor.C
    parLagrangianRedistributorFields.C

\*---------------------------------------------------------------------------*/

#ifndef parLagrangianRedistributor_H
#define parLagrangianRedistributor_H

#include "PtrList.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class mapDistributePolyMesh;
class mapDistributeBase;
class IOobjectList;
class passivePositionParticleCloud;

/*---------------------------------------------------------------------------*\
                    Class parLagrangianRedistributor Declaration
\*---------------------------------------------------------------------------*/

class parLagrangianRedistributor
{
    // Private data

        //- Source mesh reference
        const fvMesh& srcMesh_;

        //- Destination mesh reference
        const fvMesh& tgtMesh_;

        //- Distribution map reference
        const mapDistributePolyMesh& distMap_;

        //- For every src cell the target processor
        labelList destinationProcID_;

        //- For every src cell the target cell
        labelList destinationCell_;


    // Private Member Functions

        //- No copy construct
        parLagrangianRedistributor(const parLagrangianRedistributor&) = delete;

        //- No copy assignment
        void operator=(const parLagrangianRedistributor&) = delete;


public:

    // Constructors

        //- Construct from components
        parLagrangianRedistributor
        (
            const fvMesh& srcMesh,
            const fvMesh& tgtMesh,
            const label nOldCells,
            const mapDistributePolyMesh& distMap
        );


    // Member Functions

        //- Find all clouds (on all processors) and for each cloud all
        //  the objects. Result will be synchronised on all processors
        static void findClouds
        (
            const fvMesh&,
            wordList& cloudNames,
            List<wordList>& objectNames
        );

        //- Redistribute and write lagrangian positions
        autoPtr<mapDistributeBase> redistributeLagrangianPositions
        (
            passivePositionParticleCloud& cloud
        ) const;

        //- Read, redistribute and write lagrangian positions
        autoPtr<mapDistributeBase> redistributeLagrangianPositions
        (
            const word& cloudName
        ) const;

        //- Pick up any fields of a given type
        template<class Type>
        static wordList filterObjects
        (
            const IOobjectList& objects,
            const wordRes& selectedFields = wordRes()
        );

        //- Read, redistribute and write all/selected lagrangian fields
        template<class Type>
        label redistributeFields
        (
            const mapDistributeBase& map,
            const word& cloudName,
            const IOobjectList& objects,
            const wordRes& selectedFields = wordRes()
        ) const;

        //- Read, redistribute and write all/selected lagrangian fieldFields
        template<class Type>
        label redistributeFieldFields
        (
            const mapDistributeBase& map,
            const word& cloudName,
            const IOobjectList& objects,
            const wordRes& selectedFields = wordRes()
        ) const;

        //- Read and store all fields of a cloud
        template<class Container>
        static label readFields
        (
            const passivePositionParticleCloud& cloud,
            const IOobjectList& objects,
            const wordRes& selectedFields = wordRes()
        );

        //- Redistribute and write stored lagrangian fields
        template<class Container>
        label redistributeStoredFields
        (
            const mapDistributeBase& map,
            passivePositionParticleCloud& cloud
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "parLagrangianRedistributorFields.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
