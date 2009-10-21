/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
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
    Foam::cvControls

Description
    Controls for the conformalVoronoiMesh mesh generator

SourceFiles
    cvControls.C

\*---------------------------------------------------------------------------*/

#ifndef cvControls_H
#define cvControls_H

#include "IOdictionary.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class conformalVoronoiMesh;

/*---------------------------------------------------------------------------*\
                          Class cvControls Declaration
\*---------------------------------------------------------------------------*/

class cvControls
{
    // Private data

        //- Reference to the conformalVoronoiMesh holding this cvControls object
        const conformalVoronoiMesh& cvMesh_;

        //- Reference to the cvMeshDict
        const IOdictionary& cvMeshDict_;

        // General parameters

            //- Span of the domain
            scalar span_;

            //- Square of span_
            scalar spanSqr_;

            //- Activate/deactivate detailed timeChecks during meshing
            bool timeChecks_;

        // Surface conformation controls

            //- Point pair spacing coefficient - fraction of the local target
            //  cell size
            scalar pointPairDistanceCoeff_;

            //- Distance along feature edge to place mixed feature point
            //  conformation groups - to be multiple of pointPairDistanceCoeff_.
            scalar mixedFeaturePointPPDistanceCoeff_;

            //- Distance to a feature point within which surface and edge
            //  conformation points are excluded - fraction of the local target
            //  cell size
            scalar featurePointExclusionDistanceCoeff_;

            //- Distance to an existing feature edge conformation location
            //  within which other edge conformation location are excluded -
            //  fraction of the local target cell size
            scalar featureEdgeExclusionDistanceCoeff_;

            //- Surface search distance coefficient - faction of the local
            //  target cell size
            scalar surfaceSearchDistanceCoeff_;

            //- Maximum allowable protrusion through the surface before
            //  conformation points are added - fraction of the local target
            //  cell size
            scalar maxSurfaceProtrusionCoeff_;

            //- Maximum quadrant angle allowed at a concave edge before
            //  additional "mitering" lines are added
            scalar maxQuadAngle_;

            //- Now often to rebuild the surface conformation
            label surfaceConformationRebuildFrequency_;


    // Motion control controls

        //- Switch to control the output of obj files for debug
        Switch objOutput_;

        //- Number of "spokes" to use to search for secondary alignment
        //  direction
        label alignmentSearchSpokes_;

        //- cosine of angle of alignment with required direction within which a
        //  face will be accepted for rotation
        scalar cosAlignmentAcceptanceAngle_;

        // Point insertion criteria

            //- Length between Delaunay vertices above which a new Dv should be
            //  inserted - fraction of the local target cell size
            scalar insertionDistCoeff_;

            //- Minimum dual face area corresponding to long Delaunay edge where
            //  a new Dv is to be inserted - fraction of the local target cell
            //  size
            scalar faceAreaRatioCoeff_;

            //- cosine of the angle of alignment with the required direction
            //  within which a long De will be accepted for a new Dv insertion
            scalar cosInsertionAcceptanceAngle_;


        // Point removal criteria

            //- Length between Delaunay vertices below which a Dv should be
            //  removed - fraction of the local target cell size
            scalar removalDistCoeff_;


    // polyMesh filtering controls

        //- Minimum edge length allowed in the polyMesh, anything shorter will
        //  be filtered - fraction of the local target cell size
        scalar minimumEdgeLengthCoeff_;


public:

   // Constructors

        //- Construct from references to conformalVoronoiMesh and IOdictionary
        cvControls
        (
            const conformalVoronoiMesh& cvMesh,
            const IOdictionary& cvMeshDict
        );


    //- Destructor
    ~cvControls();


    // Member Functions

        // Access

            //- Return the cvMeshDict
            inline const IOdictionary& cvMeshDict() const;

            //- Return the span
            inline scalar span() const;

            //- Return spanSqr
            inline scalar spanSqr() const;

            //- Return timeChecks
            inline bool timeChecks() const;

            //- Return the pointPairDistanceCoeff
            inline scalar pointPairDistanceCoeff() const;

            //- Return the mixedFeaturePointPPDistanceCoeff
            inline scalar mixedFeaturePointPPDistanceCoeff() const;

            //- Return the featurePointExclusionDistanceCoeff
            inline scalar featurePointExclusionDistanceCoeff() const;

            //- Return the featureEdgeExclusionDistanceCoeff
            inline scalar featureEdgeExclusionDistanceCoeff() const;

            //- Return the surfaceSearchDistanceCoeff
            inline scalar surfaceSearchDistanceCoeff() const;

            //- Return the maxSurfaceProtrusionCoeff
            inline scalar maxSurfaceProtrusionCoeff() const;

            //- Return the maxQuadAngle
            inline scalar maxQuadAngle() const;

            //- Return the surfaceConformationRebuildFrequency
            inline label surfaceConformationRebuildFrequency() const;

            //- Return the objOutput Switch
            inline bool objOutput() const;

            //- Return the number of alignmentSearchSpokes to use
            inline label alignmentSearchSpokes() const;

            //- Return the cosAlignmentAcceptanceAngle
            inline scalar cosAlignmentAcceptanceAngle() const;

            //- Return the insertionDistCoeff
            inline scalar insertionDistCoeff() const;

            //- Return the faceAreaRatioCoeff
            inline scalar faceAreaRatioCoeff() const;

            //- Return the cosInsertionAcceptanceAngle
            inline scalar cosInsertionAcceptanceAngle() const;

            //- Return removalDistCoeff
            inline scalar removalDistCoeff() const;

            //- Return the minimumEdgeLengthCoeff
            inline scalar minimumEdgeLengthCoeff() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "cvControlsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
