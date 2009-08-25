/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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
    Foam::streamLineCloud

Description
    A Cloud of streamLine particles

SourceFiles
    streamLineCloud.C

\*---------------------------------------------------------------------------*/

#ifndef streamLineCloud_H
#define streamLineCloud_H

#include "Cloud.H"
#include "streamLineParticle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class streamLineCloud Declaration
\*---------------------------------------------------------------------------*/

class streamLineParticleCloud
:
    public Cloud<streamLineParticle>
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        streamLineParticleCloud(const streamLineParticleCloud&);

        //- Disallow default bitwise assignment
        void operator=(const streamLineParticleCloud&);


public:

    //- Type of parcel the cloud was instantiated for
    typedef streamLineParticle parcelType;

    // Constructors

        //- Construct given mesh
        streamLineParticleCloud
        (
            const polyMesh&,
            const word& cloudName = "defaultCloud"
        );

        //- Construct from mesh, cloud name, and a list of particles
        streamLineParticleCloud
        (
            const polyMesh& mesh,
            const word& cloudName,
            const IDLList<streamLineParticle>& particles
        );


    // Member Functions

        //- Read fields
        virtual void readFields();

        //- Write fields
        virtual void writeFields() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
