/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
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
    Foam::uniformDistance

Description

SourceFiles
    uniformDistance.C

\*---------------------------------------------------------------------------*/

#ifndef uniformDistance_H
#define uniformDistance_H

#include "cellSizeFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class uniformDistance Declaration
\*---------------------------------------------------------------------------*/

class uniformDistance
:
    public cellSizeFunction
{

private:

    // Private data

        //- Cell size
        scalar cellSize_;

        //- Distance
        scalar distance_;


public:

    //- Runtime type information
    TypeName("uniformDistance");

    // Constructors

        //- Construct from components
        uniformDistance
        (
            const dictionary& initialPointsDict,
            const conformalVoronoiMesh& cvMesh,
            const searchableSurface& surface
        );


    //- Destructor
    virtual ~uniformDistance()
    {}


    // Member Functions

        //- Modify scalar argument to the cell size specified by function.
        //  Return a boolean specifying if the function was used, i.e. false if
        //  the point was not in range of the surface for a spatially varying
        //  size.
        virtual bool cellSize(const point& pt, scalar& size) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
