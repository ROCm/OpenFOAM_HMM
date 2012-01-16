/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::tetOverlapVolume

Description
    Calculates overlap volume of two tets.

SourceFiles
    tetOverlapVolume.C

\*---------------------------------------------------------------------------*/

#ifndef tetOverlapVolume_H
#define tetOverlapVolume_H

#include "FixedList.H"
#include "labelList.H"
#include "treeBoundBox.H"

namespace Foam
{

class primitiveMesh;
class polyMesh;
class tetPoints;

/*---------------------------------------------------------------------------*\
                         Class tetOverlapVolume Declaration
\*---------------------------------------------------------------------------*/

class tetOverlapVolume
{
    // Private member functions

        //- Tet overlap
        void tetTetOverlap
        (
            const tetPoints& tetA,
            const tetPoints& tetB,
            FixedList<tetPoints, 200>& insideTets,
            label& nInside,
            FixedList<tetPoints, 200>& outsideTets,
            label& nOutside
        ) const;

        //- Tet Overlap Vol
        scalar tetTetOverlapVol
        (
            const tetPoints& tetA,
            const tetPoints& tetB
        ) const;

        //- Return a const treeBoundBox
        treeBoundBox pyrBb
        (
            const pointField& points,
            const face& f,
            const point& fc
        ) const;

public:


    //- Runtime type information
    ClassName("tetOverlapVolume");


    // Constructors

        //- Null constructor
        tetOverlapVolume();


    // Public members

        //- Return a list of cells in meshA which overlaps with cellBI in
        // meshB
        labelList overlappingCells
        (
            const polyMesh& meshA,
            const polyMesh& meshB,
            const label cellBI
        ) const;


        //- Calculates the overlap volume
        scalar cellCellOverlapVolumeMinDecomp
        (
            const primitiveMesh& meshA,
            const label cellAI,

            const primitiveMesh& meshB,
            const label cellBI,
            const treeBoundBox& cellBbB
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
