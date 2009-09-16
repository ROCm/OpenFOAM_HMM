/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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
    Foam::block

Description
    Creates a single block of cells from point coordinates, numbers of
    cells in each direction and an expansion ratio.

SourceFiles
    block.C
    blockCreate.C

\*---------------------------------------------------------------------------*/

#ifndef block_H
#define block_H

#include "pointField.H"
#include "labelList.H"

#include "blockDescriptor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Istream;
class Ostream;

/*---------------------------------------------------------------------------*\
                           Class block Declaration
\*---------------------------------------------------------------------------*/

class block
{
    // Private data

        //- Block definition
        blockDescriptor blockDef_;

        //- List of vertices
        pointField vertices_;

        //- List of cells
        labelListList cells_;

        //- Boundary patches
        labelListListList boundaryPatches_;


    // Private Member Functions

        //- Vertex label offset for a particular i,j,k position
        label vtxLabel(label i, label j, label k) const;

        //- Creates vertices for cells filling the block.
        void createPoints();

        //- Creates cells for the block.
        void createCells();

        //- Creates boundary patches for the block
        void createBoundary();

        //- Creates vertices, cells, boundary patches for the block
        void createPrimitives();

        //- Disallow default bitwise copy construct
        block(const block&);

        //- Disallow default bitwise assignment
        void operator=(const block&);

public:

    // Constructors

        //- Construct from components with Istream
        block
        (
            const pointField& blockMeshPoints,
            const curvedEdgeList&,
            Istream&
        );

        //- Construct from a block definition
        block(const blockDescriptor&);

        //- Clone
        autoPtr<block> clone() const
        {
            notImplemented("block::clone()");
            return autoPtr<block>(NULL);
        }


    // Member Functions

        // Access

            //- Return the block definition
            inline const blockDescriptor& blockDef() const
            {
                return blockDef_;
            }

            const pointField& points() const;

            const labelListList& cells() const;

            const labelListListList& boundaryPatches() const;


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const block&);
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
