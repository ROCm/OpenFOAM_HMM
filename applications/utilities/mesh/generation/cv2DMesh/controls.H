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
    Foam::controls

Description
    Controls for the 2D cv mesh generator.

SourceFiles
    controls.C

\*---------------------------------------------------------------------------*/

#ifndef controls_H
#define controls_H

#include "Switch.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class controls Declaration
\*---------------------------------------------------------------------------*/

class controls
{
    // Private data

        //- Description of data_
        const dictionary& dict_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        controls(const controls&);

        //- Disallow default bitwise assignment
        void operator=(const controls&);


public:

        //- Minimum cell size below which protusions through the surface are
        //  not split
        scalar minCellSize;

        //- Square of minCellSize
        scalar minCellSize2;

        //- The feature angle used to select corners to be
        //  explicitly represented in the mesh.
        //  0 = all features, 180 = no features
        scalar featAngle;

        //- Maximum quadrant angle allowed at a concave corner before
        //  additional "mitering" lines are added
        scalar maxQuadAngle;

        //- Should the mesh be square-dominated or of unbiased hexagons
        Switch squares;

        //- Near-wall region where cells are aligned with the wall
        scalar nearWallAlignedDist;

        //- Square of nearWallAlignedDist
        scalar nearWallAlignedDist2;

        //- Chose if the cell orientation should relax during the iterations
        //  or remain fixed to the x-y directions
        Switch relaxOrientation;

        //- Insert near-boundary point mirror or point-pairs
        Switch insertSurfaceNearestPointPairs;

        //- Mirror near-boundary points rather than insert point-pairs
        Switch mirrorPoints;

        //- Insert point-pairs vor dual-cell vertices very near the surface
        Switch insertSurfaceNearPointPairs;

        Switch writeInitialTriangulation;
        Switch writeFeatureTriangulation;
        Switch writeNearestTriangulation;
        Switch writeInsertedPointPairs;
        Switch writeFinalTriangulation;

        Switch randomiseInitialGrid;
        scalar randomPurturbation;

        label maxBoundaryConformingIter;

        //- Relaxation factor at the start of the iteration
        scalar relaxationFactorStart;

        //- Relaxation factor at the end of the iteration
        scalar relaxationFactorEnd;


    // Constructors

        controls(const dictionary& controlDict);


    //- Destructor
    ~controls();


    // Member Functions

        // Access

        // Check

        // Edit

        // Write


    // Friend Functions

    // Friend Operators

    // IOstream Operators
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#include "controlsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* // 
