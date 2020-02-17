/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Foam::unmappedPassivePositionParticleCloud

Description
    passivePositionParticleCloud but with autoMap and writing disabled.
    Only used for its objectRegistry to store lagrangian fields

SourceFiles
    unmappedPassivePositionParticleCloud.C

\*---------------------------------------------------------------------------*/

#ifndef unmappedPassivePositionParticleCloud_H
#define unmappedPassivePositionParticleCloud_H

#include "passivePositionParticleCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class unmappedPassiveParticleCloud Declaration
\*---------------------------------------------------------------------------*/

class unmappedPassivePositionParticleCloud
:
    public passivePositionParticleCloud
{

public:

    // Constructors

        //- Construct given mesh
        explicit unmappedPassivePositionParticleCloud
        (
            const polyMesh& mesh,
            const word& cloudName = cloud::defaultName,
            bool readFields = true
        )
        :
            passivePositionParticleCloud(mesh, cloudName, readFields)
        {}

        //- Construct from mesh, cloud name, and a list of particles
        unmappedPassivePositionParticleCloud
        (
            const polyMesh& mesh,
            const word& cloudName,
            const IDLList<passivePositionParticle>& particles
        )
        :
            passivePositionParticleCloud(mesh, cloudName, particles)
        {}


    //- Destructor
    virtual ~unmappedPassivePositionParticleCloud() = default;


    // Member Functions

        //- Switch off remapping of cells of particles when
        //  mesh topology changes
        virtual void autoMap(const mapPolyMesh&)
        {}

        //- Disable writing objects
        virtual bool writeObject(IOstreamOption, const bool valid) const
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
