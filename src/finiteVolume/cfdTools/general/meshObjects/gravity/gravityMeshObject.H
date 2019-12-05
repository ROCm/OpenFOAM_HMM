/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2019 OpenCFD Ltd.
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
    Foam::meshObjects::gravity

Description
    Gravitational acceleration vector
    Although termed a \em MeshObject it is registered on Time only
    and thus identical for all regions.

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
    public MeshObject
    <
        Time,
        TopologicalMeshObject,
        gravity
    >,
    public uniformDimensionedVectorField
{

public:

    //- Run-time type information
    TypeNameNoDebug("g");


    //- Construct on Time
    explicit gravity(const Time& runTime)
    :
        MeshObject<Time, TopologicalMeshObject, gravity>(runTime),
        uniformDimensionedVectorField
        (
            IOobject
            (
                "g",   // Must be identical to typeName!
                runTime.constant(),
                runTime,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE,
                false  // let MeshObject register it
            )
        )
    {}


    //- Construct on Time
    static const gravity& New(const Time& runTime)
    {
        return MeshObject<Time, TopologicalMeshObject, gravity>::New(runTime);
    }


    //- Destructor
    virtual ~gravity() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace meshObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
