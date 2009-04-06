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
    Foam::cvSearchableSurfaces

Description
    The surface geometry to be meshed with a conformalVoronoiMesh, adding
    inside/outside queries and features of the surface

SourceFiles
    cvSearchableSurfaces.C

\*---------------------------------------------------------------------------*/

#ifndef cvSearchableSurfaces_H
#define cvSearchableSurfaces_H

#include "searchableSurfaces.H"
#include "surfaceFeatures.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class conformalVoronoiMesh;

/*---------------------------------------------------------------------------*\
                     Class cvSearchableSurfaces Declaration
\*---------------------------------------------------------------------------*/

class cvSearchableSurfaces
:
    public searchableSurfaces
{
    // Private data

        //- Reference to the conformalVoronoiMesh holding this cvControls object
        const conformalVoronoiMesh& cvMesh_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        cvSearchableSurfaces(const cvSearchableSurfaces&);

        //- Disallow default bitwise assignment
        void operator=(const cvSearchableSurfaces&);


public:

    // Constructors

        //- Construct from components
        cvSearchableSurfaces
        (
            const conformalVoronoiMesh& cvMesh,
            const dictionary& geometryDict
        );

    //- Destructor
    ~cvSearchableSurfaces();


    // Member Functions

        // Access

        // Query

            //- Check if point is inside surfacesToConformTo
            bool inside(const point& pt) const;

            //- Check if point is outside surfacesToConformTo
            bool outside(const point& pt) const;

            //- Check if point is inside surfacesToConformTo by at least dist2
            bool wellInside(const point& pt, const scalar dist2) const;

            //- Check if point is outside surfacesToConformTo by at least dist2
            bool wellOutside(const point& pt, const scalar dist2) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cvSearchableSurfacesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
