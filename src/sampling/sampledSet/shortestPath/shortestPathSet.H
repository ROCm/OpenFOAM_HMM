/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::shortestPathSet

Description
    Finds shortest path (in terms of cell centres) to walk on mesh from
    any point in insidePoints to any point in outsidePoints.

Usage
    Example of function object specification:
    \verbatim
    leakFind
    {
        type            sets;

        writeControl    timeStep;
        interpolationScheme cell;
        setFormat       vtk;

        sets
        (
            leakFind
            {
                type    shortestPath;
                insidePoints   ((0.08 -0.020  -0.005) (-0.05 -0.020  -0.005));
                outsidePoints  ((-0.08 -0.020  -0.005)(0.05 -0.020  -0.005));
                axis    xyz;
            }
        );

        // Needs at least one field
        fields          ( p );
    }
    \endverbatim

SourceFiles
    shortestPathSet.C

\*---------------------------------------------------------------------------*/

#ifndef shortestPathSet_H
#define shortestPathSet_H

#include "sampledSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class topoDistanceData;

/*---------------------------------------------------------------------------*\
                       Class shortestPathSet Declaration
\*---------------------------------------------------------------------------*/

class shortestPathSet
:
    public sampledSet
{
    // Private data

        //- Originating set of points
        const pointField insidePoints_;

        //- Destination set of points
        const pointField outsidePoints_;


    // Private Member Functions

        //- Get face with least distance along route
        static label findMinFace
        (
            const polyMesh& mesh,
            const label cellI,
            const List<topoDistanceData>& allFaceInfo,
            const point& origin
        );

        //- Generate whole path
        void genSamples(const polyMesh& mesh);


public:

    //- Runtime type information
    TypeName("shortestPath");


    // Constructors

        //- Construct from components
        shortestPathSet
        (
            const word& name,
            const polyMesh& mesh,
            const meshSearch& searchEngine,
            const word& axis,
            const pointField& insidePoints,
            const pointField& outsidePoints
        );

        //- Construct from dictionary
        shortestPathSet
        (
            const word& name,
            const polyMesh& mesh,
            const meshSearch& searchEngine,
            const dictionary& dict
        );


    //- Destructor
    virtual ~shortestPathSet();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
