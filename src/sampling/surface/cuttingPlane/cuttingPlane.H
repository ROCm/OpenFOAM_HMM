/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::cuttingPlane

Description
    Constructs cutting plane through a mesh.

    No attempt at resolving degenerate cases.
    Since the cut faces can be quite ugly, they will often be triangulated.

Note
    When the cutting plane coincides with a mesh face, the cell edge on the
    positive side of the plane is taken.

SourceFiles
    cuttingPlane.C

\*---------------------------------------------------------------------------*/

#ifndef cuttingPlane_H
#define cuttingPlane_H

#include "plane.H"
#include "pointField.H"
#include "faceList.H"
#include "bitSet.H"
#include "MeshedSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class primitiveMesh;

/*---------------------------------------------------------------------------*\
                        Class cuttingPlane Declaration
\*---------------------------------------------------------------------------*/

class cuttingPlane
:
    public plane,
    public MeshedSurface<face>
{
    //- Private typedef for convenience
    typedef MeshedSurface<face> MeshStorage;


    // Private data

        //- List of cells cut by the plane
        labelList meshCells_;


    // Private Member Functions

        //- Classify sides of plane (0=BACK, 1=ONPLANE, 2=FRONT) for each point
        static PackedList<2> classifySides
        (
            const plane& pln,
            const pointField& pts
        );

        //- Determine cut cells, possibly restricted to a list of cells
        //  \return number of faces cut
        label calcCellCuts
        (
            const primitiveMesh& mesh,
            const PackedList<2>& planeSides,
            bitSet& cellCuts
        );

        //- Walk the cell cuts to create faces
        void walkCellCuts
        (
            const primitiveMesh& mesh,
            const PackedList<2>& planeSides,
            const bitSet& cellCuts,
            const bool triangulate,
            const label nFaceCuts = 0
        );


protected:

    // Constructors

        //- Construct from a plane description without any cutting
        cuttingPlane(const plane& pln);


    // Protected Member Functions

        //- Cut mesh with existing plane, restricted to a list of cells
        void performCut
        (
            const primitiveMesh& mesh,
            const bool triangulate,
            const bitSet& cellSelectionMask = bitSet()
        );

        //- Cut mesh with existing plane, restricted to a list of cells
        //  Reclaim memory for cellSelectionMask
        void performCut
        (
            const primitiveMesh& mesh,
            const bool triangulate,
            bitSet&& cellSelectionMask
        );

        //- Cut mesh with existing plane, restricted to a list of cells
        void performCut
        (
            const primitiveMesh& mesh,
            const bool triangulate,
            const labelUList& cellIdLabels
        );

        //- Remap action on triangulation or cleanup
        virtual void remapFaces(const labelUList& faceMap);


public:

    //- Debug information
    static int debug;


    // Constructors

        //- Construct from plane and mesh reference,
        //- possibly restricted to a list of cells
        cuttingPlane
        (
            const plane& pln,
            const primitiveMesh& mesh,
            const bool triangulate,
            const bitSet& cellIdLabels = bitSet()
        );

        //- Construct from plane and mesh reference,
        //- possibly restricted to a list of cells
        cuttingPlane
        (
            const plane& pln,
            const primitiveMesh& mesh,
            const bool triangulate,
            bitSet&& cellIdLabels
        );

        //- Construct from plane and mesh reference,
        //- possibly restricted to a list of cells
        cuttingPlane
        (
            const plane& pln,
            const primitiveMesh& mesh,
            const bool triangulate,
            const labelUList& cellIdLabels
        );


    // Member Functions

        //- Return the plane used
        const plane& planeDesc() const
        {
            return *this;
        }

        //- The mesh cells cut by the plane
        const labelList& meshCells() const
        {
            return meshCells_;
        }

        //- The mesh cells cut by the plane
        labelList& meshCells()
        {
            return meshCells_;
        }

        //- Have any cells been cut?
        bool cut() const
        {
            return meshCells_.size();
        }


    // Member Operators

        //- Copy assignment
        void operator=(const cuttingPlane& rhs);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
