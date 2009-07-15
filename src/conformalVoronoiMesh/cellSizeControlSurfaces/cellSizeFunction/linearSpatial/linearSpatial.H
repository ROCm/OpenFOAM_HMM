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
    Foam::linearSpatial

Description

SourceFiles
    linearSpatial.C

\*---------------------------------------------------------------------------*/

#ifndef linearSpatial_H
#define linearSpatial_H

#include "cellSizeFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class linearSpatial Declaration
\*---------------------------------------------------------------------------*/

class linearSpatial
:
    public cellSizeFunction
{

private:

    // Private data

        //- Reference point for spatial size grading
        point referencePoint_;

        //- Cell size at reference point
        scalar referenceCellSize_;

        //- Direction of cell size grading, stored as unit vector, may be
        //  supplied with any magnitude
        vector direction_;

        //- Gradient of cell size change in direction of direction_
        scalar cellSizeGradient_;


    // Private Member Functions

        //- Calculate the cell size as a function of the given position
        scalar sizeFunction(const point& pt) const;


public:

    //- Runtime type information
    TypeName("linearSpatial");

    // Constructors

        //- Construct from components
        linearSpatial
        (
            const dictionary& initialPointsDict,
            const conformalVoronoiMesh& cvMesh,
            const searchableSurface& surface
        );


    //- Destructor
    virtual ~linearSpatial()
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
