/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::meshedSurf

Description
    Abstract definition of a meshed surface defined by faces and points.

\*---------------------------------------------------------------------------*/

#ifndef meshedSurf_H
#define meshedSurf_H

#include "pointField.H"
#include "faceList.H"
#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class meshedSurf Declaration
\*---------------------------------------------------------------------------*/

class meshedSurf
{
public:

    // Constructors

        //- Construct null
        meshedSurf() = default;


    //- Destructor
    virtual ~meshedSurf() = default;


    // Access Member Functions

        //- Const access to (global) points used for the surface
        virtual const pointField& points() const = 0;

        //- Const access to the surface faces
        virtual const faceList& faces() const = 0;

        //- Const access to per-face zone/region information
        virtual const labelList& zoneIds() const = 0;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
