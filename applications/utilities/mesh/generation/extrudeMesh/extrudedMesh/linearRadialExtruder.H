/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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
    Foam::linearRadialExtruder

Description

\*---------------------------------------------------------------------------*/

#ifndef linearRadialExtruder_H
#define linearRadialExtruder_H

#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class linearRadialExtruder Declaration
\*---------------------------------------------------------------------------*/

class linearRadialExtruder
{
    // Private data

        scalar R_;


public:

    // Constructors

        //- Construct from components
        linearRadialExtruder(const scalar R)
        :
            R_(R)
        {}


    // Member Operators

        point operator()
        (
            const point& surfacePoint,
            const vector& surfaceNormal,
            const label nLayers,
            const label layer
        ) const
        {
            // radius of the surface
            scalar rs = mag(surfacePoint);
            vector rsHat = surfacePoint/rs;

            scalar delta = (R_ - rs)/nLayers;
            scalar r = rs + layer*delta;
            return r*rsHat;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
