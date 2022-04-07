/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::IOmapDistributePolyMesh

Description
    IOmapDistributePolyMesh is derived from mapDistributePolyMesh and
    IOobject to give the mapDistributePolyMesh
    automatic IO functionality via the objectRegistry.

SourceFiles
    IOmapDistributePolyMesh.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_IOmapDistributePolyMesh_H
#define Foam_IOmapDistributePolyMesh_H

#include "mapDistributePolyMesh.H"
#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class IOmapDistributePolyMesh Declaration
\*---------------------------------------------------------------------------*/

class IOmapDistributePolyMesh
:
    public regIOobject,
    public mapDistributePolyMesh
{
    // Private Member Functions

        //- Read if IOobject flags set. Return true if read.
        bool readContents();

public:

    //- Runtime type information
    TypeName("mapDistributePolyMesh");


    // Constructors

        //- Construct from IOobject
        explicit IOmapDistributePolyMesh(const IOobject& io);

        //- Construct from IOobject, copying mapDistributePolyMesh contents
        IOmapDistributePolyMesh
        (
            const IOobject& io,
            const mapDistributePolyMesh& map
        );

        //- Construct from IOobject, moving mapDistributePolyMesh contents
        IOmapDistributePolyMesh
        (
            const IOobject& io,
            mapDistributePolyMesh&& map
        );


    //- Destructor
    virtual ~IOmapDistributePolyMesh() = default;


    // Member Functions

        //- The readData method for regIOobject read operation
        virtual bool readData(Istream& is);

        //- The writeData method for regIOobject write operation
        virtual bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
