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
    Foam::surfaceOffsetLinearDistance

Description

SourceFiles
    surfaceOffsetLinearDistance.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceOffsetLinearDistance_H
#define surfaceOffsetLinearDistance_H

#include "cellSizeFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class surfaceOffsetLinearDistance Declaration
\*---------------------------------------------------------------------------*/

class surfaceOffsetLinearDistance
:
    public cellSizeFunction
{

private:

    // Private data

        //- cell size at the surface
        scalar surfaceCellSize_;

        //- cell size at distance_ from the surface
        scalar distanceCellSize_;

        //- Offset distance from surface for constant size portion
        scalar surfaceOffset_;

        //- Total distance from the surface to control over (distance +
        //  surfaceOffset)
        scalar totalDistance_;

        //- totalDistance squared
        scalar totalDistanceSqr_;

        //- storing gradient for linear function
        scalar gradient_;

        //- storing intercept for linear function
        scalar intercept_;


    // Private Member Functions

        //- Calculate the cell size as a function of the given distance
        scalar sizeFunction(scalar d) const;


public:

    //- Runtime type information
    TypeName("surfaceOffsetLinearDistance");

    // Constructors

        //- Construct from components
        surfaceOffsetLinearDistance
        (
            const dictionary& initialPointsDict,
            const conformalVoronoiMesh& cvMesh,
            const searchableSurface& surface
        );


    //- Destructor
    virtual ~surfaceOffsetLinearDistance()
    {}


    // Member Functions

        //- Modify scalar argument to the cell size specified by function.
        //  Return a boolean specifying if the function was used, i.e. false if
        //  the point was not in range of the surface for a spatially varying
        //  size.
        virtual bool cellSize
        (
            const point& pt,
            scalar& size,
            bool isSurfacePoint = false
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
