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
    Constructs plane through mesh.

    No attempt at resolving degenerate cases. Since the cut faces are
    usually quite ugly, they will always be triangulated.

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
#include "MeshedSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

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

        //- Determine cut cells, possibly restricted to a list of cells
        void calcCutCells
        (
            const primitiveMesh&,
            const scalarField& dotProducts,
            const labelUList& cellIdLabels = labelUList::null()
        );

        //- Determine intersection points (cutPoints).
        void intersectEdges
        (
            const primitiveMesh&,
            const scalarField& dotProducts,
            List<label>& edgePoint
        );

        //- Walk circumference of cell, starting from startEdgeI crossing
        //  only cut edges. Record cutPoint labels in faceVerts.
        static bool walkCell
        (
            const primitiveMesh&,
            const labelUList& edgePoint,
            const label celli,
            const label startEdgei,
            DynamicList<label>& faceVerts
        );

        //- Determine cuts for all cut cells.
        void walkCellCuts
        (
            const primitiveMesh& mesh,
            const bool triangulate,
            const labelUList& edgePoint
        );


protected:

    // Constructors

        //- Construct plane description without cutting
        cuttingPlane(const plane& pln);


    // Protected Member Functions

        //- Recut mesh with existing plane, restricted to a list of cells
        void reCut
        (
            const primitiveMesh&,
            const bool triangulate,
            const labelUList& cellIdLabels = labelUList::null()
        );

        //- Remap action on triangulation or cleanup
        virtual void remapFaces(const labelUList& faceMap);


public:

    // Constructors

        //- Construct from plane and mesh reference,
        //- possibly restricted to a list of cells
        cuttingPlane
        (
            const plane& pln,
            const primitiveMesh& mesh,
            const bool triangulate,
            const labelUList& cellIdLabels = labelUList::null()
        );


    // Member Functions

        //- Return the plane used
        const plane& planeDesc() const
        {
            return static_cast<const plane&>(*this);
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
