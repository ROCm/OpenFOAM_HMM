/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2011 OpenCFD Ltd.
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
    Foam::parallelHierarchicalDensityWeightedStochastic

Description
    Choose random points inside the domain and place them with a probability
    proportional to the target density of points.

SourceFiles
    parallelHierarchicalDensityWeightedStochastic.C

\*---------------------------------------------------------------------------*/

#ifndef parallelHierarchicalDensityWeightedStochastic_H
#define parallelHierarchicalDensityWeightedStochastic_H

#include "initialPointsMethod.H"

#include "fvMesh.H"
#include "hexRef8.H"
#include "cellSet.H"
#include "meshTools.H"
#include "polyTopoChange.H"
#include "mapPolyMesh.H"
#include "decompositionMethod.H"
#include "fvMeshDistribute.H"
#include "mapDistributePolyMesh.H"
#include "globalIndex.H"
#include "treeBoundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class parallelHierarchicalDensityWeightedStochastic Declaration
\*---------------------------------------------------------------------------*/

class parallelHierarchicalDensityWeightedStochastic
:
    public initialPointsMethod
{

private:

    // Private data

        //- Trial points attempted to be placed in all boxes
        mutable label globalTrialPoints_;

        //- Smallest minimum cell size allowed, i.e. to avoid high initial
        //  population of areas of small size
        scalar minCellSizeLimit_;

        //- Maximum normal level of recursion, can be more if a high density
        //  ratio is detected
        label minLevels_;

        //- Maximum allowed ratio of cell size in a box
        scalar maxSizeRatio_;

        //- How fine should the initial sample of the volume a box be to
        //  investigate it cell sizes and volume fraction
        label volRes_;

        //- How fine should the initial sample of the surface of a box be to
        //  investigate if it is near to a the geometry.
        label surfRes_;

    // Private Member Functions

        //- Print details of the decomposed mesh
        void printMeshData(const polyMesh& mesh) const;

        //- Estimate the number of vertices that will be in this cell, returns
        //  true if the cell is to be split because of the density ratio inside
        //  it
        bool estimateCellWeight
        (
            const polyMesh& mesh,
            label cellI,
            label volType,
            scalar& weightEstimate
        ) const;

        //- Select cells for refinement at the surface of the geometry to be
        //  meshed
        labelList selectRefinementCells
        (
            const hexRef8& meshCutter,
            labelList& volumeStatus,
            volScalarField& cellWeights
        ) const;

public:

    //- Runtime type information
    TypeName("parallelHierarchicalDensityWeightedStochastic");

    // Constructors

        //- Construct from components
        parallelHierarchicalDensityWeightedStochastic
        (
            const dictionary& initialPointsDict,
            const conformalVoronoiMesh& cvMesh
        );


    //- Destructor
    virtual ~parallelHierarchicalDensityWeightedStochastic()
    {}


    // Member Functions

        //- Return the initial points for the conformalVoronoiMesh
        virtual std::vector<Vb::Point> initialPoints() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
