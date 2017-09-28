/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::passivePositionParticleCloud

Description
    A Cloud of passive position particles

SourceFiles
    passivePositionParticleCloud.C

\*---------------------------------------------------------------------------*/

#ifndef passivePositionParticleCloud_H
#define passivePositionParticleCloud_H

#include "Cloud.H"
#include "passivePositionParticle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class passivePositionParticleCloud Declaration
\*---------------------------------------------------------------------------*/

class passivePositionParticleCloud
:
    public Cloud<passivePositionParticle>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        passivePositionParticleCloud(const passivePositionParticleCloud&);

        //- Disallow default bitwise assignment
        void operator=(const passivePositionParticleCloud&);


public:

    // Constructors

        //- Construct given mesh
        passivePositionParticleCloud
        (
            const polyMesh&,
            const word& cloudName = "defaultCloud",
            bool readFields = true
        );

        //- Construct from mesh, cloud name, and a list of particles
        passivePositionParticleCloud
        (
            const polyMesh& mesh,
            const word& cloudName,
            const IDLList<passivePositionParticle>& particles
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
