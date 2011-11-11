/*--------------------------------*- C++ -*----------------------------------*\
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
    Foam::patchTo2DpolyMesh

Description
    Convert a primitivePatch into a 2D polyMesh.

SourceFiles
    patchTo2DpolyMeshI.H
    patchTo2DpolyMesh.C
    patchTo2DpolyMeshIO.C

\*---------------------------------------------------------------------------*/

#ifndef patchTo2DpolyMesh_H
#define patchTo2DpolyMesh_H

#include "EdgeMap.H"
#include "polyMesh.H"
#include "MeshedSurface.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class patchTo2DpolyMesh Declaration
\*---------------------------------------------------------------------------*/

class patchTo2DpolyMesh
{
    // Private data

        const MeshedSurface<face>& patch_;

        const wordList& patchNames_;

        const labelList& patchSizes_;

        labelList patchStarts_;

        const EdgeMap<label>& mapEdgesRegion_;

        pointField points_;
        faceList faces_;
        labelList owner_;
        labelList neighbour_;

        //- Description of data_
        void flipFaceOrder();

        void createNeighbours();

        labelList internalFaceOrder();

        void addPatchFacesToFaces();

        void addPatchFacesToOwner();

        void createPolyMeshComponents();


    // Private Member Functions

        //- Disallow default bitwise copy construct
        patchTo2DpolyMesh(const patchTo2DpolyMesh&);

        //- Disallow default bitwise assignment
        void operator=(const patchTo2DpolyMesh&);


public:

    // Constructors

        //- Construct from a primitivePatch
        patchTo2DpolyMesh
        (
            const MeshedSurface<face>& patch,
            const wordList& patchNames,
            const labelList& patchSizes,
            const EdgeMap<label>& mapEdgesRegion
        );


    //- Destructor
    ~patchTo2DpolyMesh();


    // Member Functions

        // Access
        pointField& points()
        {
            return points_;
        }

        faceList& faces()
        {
            return faces_;
        }

        labelList& owner()
        {
            return owner_;
        }

        labelList& neighbour()
        {
            return neighbour_;
        }

        const wordList& patchNames() const
        {
            return patchNames_;
        }

        const labelList& patchSizes() const
        {
            return patchSizes_;
        }

        const labelList& patchStarts() const
        {
            return patchStarts_;
        }
        // Check

        // Edit
        void createMesh();

        // Write

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* // 
