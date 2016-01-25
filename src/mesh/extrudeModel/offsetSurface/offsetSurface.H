/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014 OpenFOAM Foundation
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
    Foam::extrudeModels::offsetSurface

Description
    Extrudes by interpolating points from one surface to the other. Surfaces
    have to be topologically identical i.e. one has to be an offsetted version
    of the other.

\*---------------------------------------------------------------------------*/

#ifndef offsetSurface_H
#define offsetSurface_H

#include "point.H"
#include "extrudeModel.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class triSurface;
class triSurfaceSearch;


namespace extrudeModels
{

/*---------------------------------------------------------------------------*\
                    Class offsetSurface Declaration
\*---------------------------------------------------------------------------*/

class offsetSurface
:
    public extrudeModel
{
    // Private data

        //- surface
        autoPtr<triSurface> baseSurfPtr_;

        //- search engine
        autoPtr<triSurfaceSearch> baseSearchPtr_;

        //- offsets
        autoPtr<triSurface> offsetSurfPtr_;

        //- search engine
        autoPtr<triSurfaceSearch> offsetSearchPtr_;

        // Whether to re-project onto offsetted surface
        const Switch project_;

public:

    //- Runtime type information
    TypeName("offsetSurface");

    // Constructors

        //- Construct from dictionary
        offsetSurface(const dictionary& dict);


    //- Destructor
    virtual ~offsetSurface();


    // Member Operators

        //- Return point
        point operator()
        (
            const point& surfacePoint,
            const vector& surfaceNormal,
            const label layer
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace extrudeModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
