/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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
    Foam::pointBoundaryMesh

Description
    Foam::pointBoundaryMesh

SourceFiles
    pointBoundaryMesh.C

\*---------------------------------------------------------------------------*/

#ifndef pointBoundaryMesh_H
#define pointBoundaryMesh_H

#include "pointPatchList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class pointMesh;
class polyBoundaryMesh;
class wordRes;

/*---------------------------------------------------------------------------*\
                       Class pointBoundaryMesh Declaration
\*---------------------------------------------------------------------------*/

class pointBoundaryMesh
:
    public pointPatchList
{
    // Private Data

        //- Reference to mesh
        const pointMesh& mesh_;


    // Private Member Functions

        //- Calculate geometry for the patches (transformation tensors etc.)
        void calcGeometry();

        //- No copy construct
        pointBoundaryMesh(const pointBoundaryMesh&) = delete;

        //- No copy assignment
        void operator=(const pointBoundaryMesh&) = delete;


public:

    //- Declare friendship with pointMesh
    friend class pointMesh;


    // Constructors

        //- Construct from polyBoundaryMesh
        pointBoundaryMesh
        (
            const pointMesh&,
            const polyBoundaryMesh&
        );


    // Member Functions

        //- Return the mesh reference
        const pointMesh& mesh() const noexcept
        {
            return mesh_;
        }

        //- Return (sorted) patch indices for all matches.
        //  A no-op (returns empty list) for an empty matcher
        labelList indices(const wordRe& matcher, const bool useGroups) const;

        //- Return (sorted) patch indices for all matches.
        //  A no-op (returns empty list) for an empty matcher
        labelList indices(const wordRes& matcher, const bool useGroups) const;

        //- Find patch index given a name
        //  A no-op (returns -1) for an empty patchName
        label findPatchID(const word& patchName) const;

        //- Correct polyBoundaryMesh after moving points
        void movePoints(const pointField&);

        //- Correct polyBoundaryMesh after topology update
        void updateMesh();


    // Housekeeping

        //- Identical to the indices() method (AUG-2018)
        FOAM_DEPRECATED_FOR(2018-08, "indices() method")
        labelList findIndices(const wordRe& key, bool useGroups) const
        {
            return indices(key, useGroups);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
