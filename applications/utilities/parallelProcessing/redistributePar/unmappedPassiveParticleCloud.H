/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd.
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
    Foam::unmappedPassiveParticleCloud

Description
    passiveParticleCloud but with autoMap and writing disabled. Only used
    for its objectRegistry to store lagrangian fields

SourceFiles
    unmappedPassiveParticleCloud.C

\*---------------------------------------------------------------------------*/

#ifndef unmappedPassiveParticleCloud_H
#define unmappedPassiveParticleCloud_H

#include "passiveParticleCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class unmappedPassiveParticleCloud Declaration
\*---------------------------------------------------------------------------*/

class unmappedPassiveParticleCloud
:
    public passiveParticleCloud
{

public:

    // Constructors

        //- Construct given mesh
        unmappedPassiveParticleCloud
        (
            const polyMesh& mesh,
            const word& cloudName = "defaultCloud",
            bool readFields = true
        )
        :
            passiveParticleCloud(mesh, cloudName, readFields)
        {}

        //- Construct from mesh, cloud name, and a list of particles
        unmappedPassiveParticleCloud
        (
            const polyMesh& mesh,
            const word& cloudName,
            const IDLList<passiveParticle>& particles
        )
        :
            passiveParticleCloud(mesh, cloudName, particles)
        {}

        //- Destructor
        virtual ~unmappedPassiveParticleCloud()
        {}

        //- Switch off remapping of cells of particles when
        //  mesh topology changes
        virtual void autoMap(const mapPolyMesh&)
        {}

        //- Switch off writing the objects
        virtual bool writeObject
        (
            IOstream::streamFormat fmt,
            IOstream::versionNumber ver,
            IOstream::compressionType cmp
        ) const
        {
            return true;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
