/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::blockEdges::arcEdge

Description
    Defines the arcEdge of a circle in terms of 3 points on its circumference

SourceFiles
    arcEdge.C

\*---------------------------------------------------------------------------*/

#ifndef blockEdges_arcEdge_H
#define blockEdges_arcEdge_H

#include "blockEdge.H"
#include "cylindricalCS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blockEdges
{

/*---------------------------------------------------------------------------*\
                           Class arcEdge Declaration
\*---------------------------------------------------------------------------*/

class arcEdge
:
    public blockEdge
{
    // Private data

        // Begin, mid, end points
        point p1_, p2_, p3_;

        //- The arc angle (in degrees)
        scalar angle_;

        //- The arc radius
        scalar radius_;

        //- The local cylindrical coordinate system
        coordSystem::cylindrical cs_;


    // Private Member Functions

        //- Calculate the angle, radius and axis
        //  \return the cylindrical coordinate system
        coordSystem::cylindrical calcAngle();

        //- No copy construct
        arcEdge(const arcEdge&) = delete;

        //- No copy assignment
        void operator=(const arcEdge&) = delete;


public:

    //- Runtime type information
    TypeName("arc");


    // Constructors

        //- Construct from components
        arcEdge
        (
            const pointField& points,
            const label start,
            const label end,
            const point& pMid
        );

        //- Construct from Istream setting pointsList
        arcEdge
        (
            const dictionary& dict,
            const label index,
            const searchableSurfaces& geometry,
            const pointField& points,
            Istream&
        );


    //- Destructor
    virtual ~arcEdge() = default;


    // Member Functions

        //- The point corresponding to the curve parameter [0-1]
        point position(const scalar lambda) const;

        //- The length of the curve
        scalar length() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace blockEdges
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
