/*---------------------------------------------------------------------------*\
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
    Foam::backgroundMeshDecomposition

Description
    Store a background polyMesh to use for the decomposition of space and
    queries for parallel conformalVoronoiMesh.

    The requirements are:

    + To have a decomposition of space which can quickly interrogate an
      arbitrary location from any processor to reliably and unambiguously
      determine which processor owns the space that the point is in, i.e. as
      the vertices move, or need inserted as part of the surface conformation,
      send them to the correct proc.

    + To be able to be dynamically built, refined and redistributed to other
      procs the partitioning as the meshing progresses to balance the load.

    + To be able to query whether a sphere (the circumsphere of a Delaunay tet)
      overlaps any part of the space defined by the structure, and whether a
      ray (Voronoi edge) penetrates any part of the space defined by the
      structure, this is what determines if points get referred to a processor.

SourceFiles
    backgroundMeshDecompositionI.H
    backgroundMeshDecomposition.C

\*---------------------------------------------------------------------------*/

#ifndef backgroundMeshDecomposition_H
#define backgroundMeshDecomposition_H

#include "conformalVoronoiMesh.H"
#include "fvMesh.H"
#include "hexRef8.H"
#include "cellSet.H"
#include "meshTools.H"
#include "polyTopoChange.H"
#include "mapPolyMesh.H"
#include "decompositionMethod.H"
#include "fvMeshDistribute.H"
#include "removeCells.H"
#include "mapDistributePolyMesh.H"
#include "globalIndex.H"
#include "treeBoundBox.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class backgroundMeshDecomposition Declaration
\*---------------------------------------------------------------------------*/

class backgroundMeshDecomposition
{
    // Private data

        //- Method details dictionary
        dictionary coeffsDict_;

        //- Reference to the conformalVoronoiMesh holding this object
        const conformalVoronoiMesh& cvMesh_;

        //- Mesh stored on for this processor, specifiying the domain that it
        //  is responsible for.
        fvMesh mesh_;

        //- Scale of a cell span vs cell size used to decide to refine a cell
        scalar spanScale_;

        //- Smallest minimum cell size allowed, i.e. to avoid high initial
        //  refinement of areas of small size
        scalar minCellSizeLimit_;

        //- Minimum normal level of refinement
        label minLevels_;

        //- How fine should the initial sample of the volume a box be to
        //  investigate the local cell size
        label volRes_;


    // Private Member Functions

        void initialRefinement();

        //- Print details of the decomposed mesh
        void printMeshData(const polyMesh& mesh) const;

        //- Estimate the number of vertices that will be in this cell, returns
        //  true if the cell is to be split because of the density ratio inside
        //  it
        bool refineCell
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

        //- Disallow default bitwise copy construct
        backgroundMeshDecomposition(const backgroundMeshDecomposition&);

        //- Disallow default bitwise assignment
        void operator=(const backgroundMeshDecomposition&);


public:

    //- Runtime type information
    ClassName("backgroundMeshDecomposition");


    // Constructors

        //- Construct from components
        backgroundMeshDecomposition
        (
            const dictionary& coeffsDict,
            const conformalVoronoiMesh& cvMesh
        );


    //- Destructor
    ~backgroundMeshDecomposition();


    // Member Functions

    void distribute() const;

    //- Check if the point is in the domain handled by this processor
    bool positionOnThisProc(const point& pt) const;

    //- Which processor's domain handles this point
    label positionProc(const point& pt) const;


        // Access

        // Check

        // Edit

        // Write

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "backgroundMeshDecompositionI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
