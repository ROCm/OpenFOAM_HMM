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

    For a dictionary specification:
    \table
        Property | Description                             | Required | Default
        type     | shortestPath                            | yes      |
        axis     | x, y, z, xyz, distance                  | yes      |
        insidePoints  | The inside points                  | yes      |
        outsidePoints | The outside points                 | yes      |
    \endtable

SourceFiles
    shortestPathSet.C

\*---------------------------------------------------------------------------*/

#ifndef shortestPathSet_H
#define shortestPathSet_H

#include "sampledSet.H"
#include "PackedBoolList.H"

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
            const PackedBoolList& isLeakPoint,
            const bool minDistance,
            const point& origin
        );

        //- Sync the leak path data
        void sync
        (
            const polyMesh& mesh,
            PackedBoolList& isLeakFace,
            PackedBoolList& isLeakPoint,
            const label celli,
            point& origin,
            bool& findMinDistance
        ) const;

        //- Calculate (topological) distance to cellI
        void calculateDistance
        (
            const label iter,
            const polyMesh& mesh,
            const label cellI,

            List<topoDistanceData>& allFaceInfo,
            List<topoDistanceData>& allCellInfo
        ) const;

        //- Checks if face uses a leak point
        bool touchesWall
        (
            const polyMesh& mesh,
            const label facei,

            PackedBoolList& isLeakFace,
            const PackedBoolList& isLeakPoint
        ) const;

        //- Calculate path between insideCelli (-1 if not on current processor)
        //  and outsidePoint. Appends cellcentres on path to track.
        //      isLeakCell  : track passes through cell
        //      isLeakFace  : faces of leak cells that are also on boundary
        //      isLeakPoint : points  of leak faces       ,,
        void genSamples
        (
            const bool markLeakPath,
            const label maxIter,
            const polyMesh& mesh,
            const boolList& isBlockedFace,
            const point& insidePoint,
            const label insideCelli,
            const point& outsidePoint,

            DynamicList<point>& samplingPts,
            DynamicList<label>& samplingCells,
            DynamicList<label>& samplingFaces,
            DynamicList<label>& samplingSegments,
            DynamicList<scalar>& samplingCurveDist,
            PackedBoolList& isLeakCell,
            PackedBoolList& isLeakFace,
            PackedBoolList& isLeakPoint
        ) const;

        //- Generate whole path. With markLeakPath it will block all faces
        //  along the whole path so will maximise the chances of finding
        //  multiple gaps. With markLeakPath=false it will only block any
        //  faces connected to a boundary. This makes for the nicest
        //  hole-filling.
        void genSamples
        (
            const bool markLeakPath,    // mark all cells along path
            const label maxIter,
            const polyMesh& mesh,
            const labelUList& wallPatches,
            const boolList& blockedFace
        );


public:

    //- Runtime type information
    TypeName("shortestPath");


    // Constructors

        //- Construct from components. blockedFace is an optional specification
        //  of face that behave as if a wall
        shortestPathSet
        (
            const word& name,
            const polyMesh& mesh,
            const meshSearch& searchEngine,
            const word& axis,
            const bool markLeakPath,
            const label maxIter,
            const labelUList& wallPatches,
            const pointField& insidePoints,
            const pointField& outsidePoints,
            const boolList& blockedFace
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
    virtual ~shortestPathSet() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
