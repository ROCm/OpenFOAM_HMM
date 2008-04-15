/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::mapDistributeLagrangian

Description
    Class containing mesh-to-mesh mapping information for particles

SourceFiles
    mapDistributeLagrangian.C

\*---------------------------------------------------------------------------*/

#ifndef mapDistributeLagrangian_H
#define mapDistributeLagrangian_H

#include "mapDistributePolyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                           Class mapDistributeLagrangian Declaration
\*---------------------------------------------------------------------------*/

class mapDistributeLagrangian
{
    // Private data

        const label nNewParticles_;

        //- Maps from subsetted mesh back to original mesh
        const labelListList subParticleMap_;

        //- Maps from subsetted mesh to new reconstructed mesh
        const labelListList constructParticleMap_;

        //- Per element in subsetted mesh the cell label
        const labelListList constructCellLabels_;


public:

    // Constructors

        //- Construct from components
        mapDistributeLagrangian
        (
            const label nNewParticles,
            const labelListList& subParticleMap,
            const labelListList& constructParticleMap,
            const labelListList& constructCellLabels
        )
        :
            nNewParticles_(nNewParticles),
            subParticleMap_(subParticleMap),
            constructParticleMap_(constructParticleMap),
            constructCellLabels_(constructCellLabels)
        {}

        //- Construct from components and steal storage
        mapDistributeLagrangian
        (
            const label nNewParticles,
            labelListList& subParticleMap,
            labelListList& constructParticleMap,
            labelListList& constructCellLabels,
            const bool reUse
        )
        :
            nNewParticles_(nNewParticles),
            subParticleMap_(subParticleMap, reUse),
            constructParticleMap_(constructParticleMap, reUse),
            constructCellLabels_(constructCellLabels, reUse)
        {}


    // Member Functions

        // Access

            //- From subsetted lagrangian particles back to original particles
            const labelListList& subParticleMap() const
            {
                return subParticleMap_;
            }

            //- From received particles to new reconstructed particles
            const labelListList& constructParticleMap() const
            {
                return constructParticleMap_;
            }

            //- Per received particle the destination cell label
            const labelListList& constructCellLabels() const
            {
                return constructCellLabels_;
            }


        // Edit

            //- distribute list of lagrangian data
            template<class T>
            void distributeLagrangianData(List<T>& lst) const
            {
                mapDistributePolyMesh::distribute
                (
                    nNewParticles_,
                    subParticleMap(),
                    constructParticleMap(),
                    lst
                );
            }

            //- Correct for topo change.
            void updateMesh(const mapPolyMesh&)
            {
                notImplemented
                (
                    "mapDistributeLagrangian::updateMesh(const mapPolyMesh&)"
                );
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
