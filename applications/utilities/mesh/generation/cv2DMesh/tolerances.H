/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2011 OpenCFD Ltd.
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
    Foam::tolerances

Description
    Tolerances for the CV 2D mesher.

SourceFiles
    tolerances.C

\*---------------------------------------------------------------------------*/

#ifndef tolerances_H
#define tolerances_H

#include "dictionary.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class tolerances Declaration
\*---------------------------------------------------------------------------*/

class tolerances
{
    // Private data

        //- Description of data_


    // Private Member Functions

        //- Disallow default bitwise copy construct
        tolerances(const tolerances&);

        //- Disallow default bitwise assignment
        void operator=(const tolerances&);


public:

        //- Maximum cartesian span of the geometry
        scalar span;

        //- Square of span
        scalar span2;

        //- Minumum edge-length of the cell size below which protusions
        //  through the surface are not split
        scalar minEdgeLen;

        //- Square of minEdgeLen
        scalar minEdgeLen2;

        //- Maximum notch size below which protusions into the surface are
        //  not filled
        scalar maxNotchLen;

        //- Square of maxNotchLen
        scalar maxNotchLen2;

        //- The minimum distance alowed between a dual-cell vertex
        // and the surface before a point-pair is introduced
        scalar minNearPointDist;

        //- Square of minNearPoint
        scalar minNearPointDist2;

        //- Distance between boundary conforming point-pairs
        scalar ppDist;

        //- Square of ppDist
        scalar ppDist2;


    // Constructors

        //- Construct null
        tolerances
        (
            const dictionary& controlDict,
            scalar minCellSize,
            const boundBox&
        );


    //- Destructor
    ~tolerances();


    // Member Functions

        // Access

        // Check

        // Edit

        // Write

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* // 
