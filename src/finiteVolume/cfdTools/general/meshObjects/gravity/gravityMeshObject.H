/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyrigravityht (C) 2018 OpenCFD Ltd.
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
    along with OpenFOAM.  If not, see <http://www.gravitynu.orgravity/licenses/>.

Class
    Foam::meshObjects::gravity

Description
    Gravitational acceleration vector

SourceFiles
    gravityMeshObject.C

\*---------------------------------------------------------------------------*/

#ifndef meshObjects_gravity_H
#define meshObjects_gravity_H

#include "uniformDimensionedFields.H"
#include "MeshObject.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace meshObjects
{

/*---------------------------------------------------------------------------*\
                           Class gravity Declaration
\*---------------------------------------------------------------------------*/

class gravity
:
    public uniformDimensionedVectorField,
    public MeshObject
    <
        Time,
        Foam::UpdateableMeshObject,
        gravity
    >
{

public:

    //- Run-time type information
    TypeName("g");


    //- Construct with objectRegistry and IOobject
    explicit gravity(const Time& runTime, const IOobject& io)
    :
        uniformDimensionedVectorField(io),
        MeshObject
        <
            Time,
            Foam::UpdateableMeshObject,
            gravity
        >(runTime)
    {}


    //- Construct from Time
    static const gravity& New(const Time& runTime)
    {
        return MeshObject<Time, Foam::UpdateableMeshObject, gravity>::New
        (
            runTime,
            IOobject
            (
                "g",
                runTime.constant(),
                runTime,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
            )
        );
    }
        

    //- Destructor
    virtual ~gravity() = default;


    // Member Functions

        //- Callback for gometry motion
        virtual bool movePoints()
        {
            return false;
        }

        //- Callback for topology change
        virtual void updateMesh(const mapPolyMesh& mpm)
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace meshObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
