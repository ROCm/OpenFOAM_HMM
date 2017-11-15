/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd.
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
    Foam::refinementParameters

Description
    Simple container to keep together refinement specific information.

SourceFiles
    refinementParameters.C

\*---------------------------------------------------------------------------*/

#ifndef refinementParameters_H
#define refinementParameters_H

#include "dictionary.H"
#include "pointField.H"
#include "Switch.H"
#include "wordPairHashTable.H"
#include "surfaceZonesInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class polyMesh;

/*---------------------------------------------------------------------------*\
                           Class refinementParameters Declaration
\*---------------------------------------------------------------------------*/

class refinementParameters
{
    // Private data

        //- Total number of cells
        const label maxGlobalCells_;

        //- Per processor max number of cells
        const label maxLocalCells_;

        //- When to stop refining
        const label minRefineCells_;

        //- Curvature
        scalar curvature_;

        //- Planarity criterion
        scalar planarAngle_;

        //- Number of layers between different refinement levels
        const label nBufferLayers_;


        // Selection of areas

            //- Areas not to keep
            const pointField locationsOutsideMesh_;

            //- Areas to keep
            pointField locationsInMesh_;

            //- Region for location
            wordList zonesInMesh_;

            //- Information on how to handle faces on faceZones
            dictionary faceZoneControls_;


        //- FaceZone faces allowed which have owner and neighbour in same
        //  cellZone?
        Switch allowFreeStandingZoneFaces_;

        //- Use old topology based problem-cell removal (cells with 8 points
        //  on surface)
        Switch useTopologicalSnapDetection_;

        //- Allowed load unbalance
        scalar maxLoadUnbalance_;

        Switch handleSnapProblems_;

        Switch interfaceRefine_;

        label nErodeCellZone_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        refinementParameters(const refinementParameters&);

        //- Disallow default bitwise assignment
        void operator=(const refinementParameters&);


public:

    // Constructors

        //- Construct from dictionary - new syntax
        refinementParameters(const dictionary& dict);


    // Member Functions

        // Access

            //- Total number of cells
            label maxGlobalCells() const
            {
                return maxGlobalCells_;
            }

            //- Per processor max number of cells
            label maxLocalCells() const
            {
                return maxLocalCells_;
            }

            //- When to stop refining
            label minRefineCells() const
            {
                return minRefineCells_;
            }

            //- Curvature
            scalar curvature() const
            {
                return curvature_;
            }

            //- Angle when two intersections are considered to be planar
            scalar planarAngle() const
            {
                return planarAngle_;
            }

            //- Number of layers between different refinement levels
            label nBufferLayers() const
            {
                return nBufferLayers_;
            }

            //- Areas to keep
            const pointField& locationsInMesh() const
            {
                return locationsInMesh_;
            }

            //- Per area the zone name
            const wordList& zonesInMesh() const
            {
                return zonesInMesh_;
            }

            //- Optional points which are checked to be outside the mesh
            const pointField& locationsOutsideMesh() const
            {
                return locationsOutsideMesh_;
            }

            //- Are zone faces allowed only inbetween different cell zones
            //  or also just free standing?
            bool allowFreeStandingZoneFaces() const
            {
                return allowFreeStandingZoneFaces_;
            }

            //- Use old topology based problem-cell removal
            //  (cells with 8 points on surface)
            bool useTopologicalSnapDetection() const
            {
                return useTopologicalSnapDetection_;
            }

            //- Allowed load unbalance
            scalar maxLoadUnbalance() const
            {
                return maxLoadUnbalance_;
            }

            bool handleSnapProblems() const
            {
                return handleSnapProblems_;
            }

            //- Refine cell with opposite faces with different refinement level
            bool interfaceRefine() const
            {
                return interfaceRefine_;
            }

            label nErodeCellZone() const
            {
                return nErodeCellZone_;
            }


        // Other

            //- Get patchInfo and faceType for faceZone
            dictionary getZoneInfo
            (
                const word& fzName,
                surfaceZonesInfo::faceZoneType& faceType
            ) const;

            //- Add cellZones to mesh. Return indices of cellZones (or -1)
            labelList addCellZonesToMesh(polyMesh&) const;

            //- Checks that cells are in mesh. Returns cells (or -1) they
            //  are in.
            static labelList findCells
            (
                const bool checkInsideMesh,
                const polyMesh&,
                const pointField& locations
            );

            //- Extract indices of named locations (so excludes 'keepPoints')
            static labelList zonedLocations(const wordList& zonesInMesh);

            //- Extract indices of unnamed locations ('keepPoints')
            static labelList unzonedLocations(const wordList& zonesInMesh);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
