/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::injectedParticleCloud

Group

Description

SourceFiles
    injectedParticleCloud.C

\*---------------------------------------------------------------------------*/

#ifndef injectedParticleCloud_H
#define injectedParticleCloud_H

#include "Cloud.H"
#include "injectedParticle.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


/*---------------------------------------------------------------------------*\
                   Class injectedParticleCloud Declaration
\*---------------------------------------------------------------------------*/

class injectedParticleCloud
:
    public Cloud<injectedParticle>
{

private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        injectedParticleCloud(const injectedParticleCloud&);

        //- Disallow default bitwise assignment
        void operator=(const injectedParticleCloud&);


public:

    // Constructors

        //- Construct from mesh and cloud name
        injectedParticleCloud
        (
            const polyMesh& mesh,
            const word& name,
            const bool readFields = true
        );

        //- Copy constructor with new name
        injectedParticleCloud(const injectedParticleCloud& c, const word& name);

        //- Construct and return clone based on (this) with new name
        virtual autoPtr<injectedParticleCloud> clone(const word& name) const
        {
            return autoPtr<injectedParticleCloud>
            (
                new injectedParticleCloud(*this, name)
            );
        }


    //- Destructor
    virtual ~injectedParticleCloud();


    // Member Functions

        // I-O

            //- Write particle fields as objects into the obr registry
            virtual void writeObjects(objectRegistry& obr) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
