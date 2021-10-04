/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::fvBoundaryMesh

Description
    Foam::fvBoundaryMesh

SourceFiles
    fvBoundaryMesh.C

\*---------------------------------------------------------------------------*/

#ifndef fvBoundaryMesh_H
#define fvBoundaryMesh_H

#include "fvPatchList.H"
#include "lduInterfacePtrsList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class fvMesh;
class polyBoundaryMesh;
class wordRes;

/*---------------------------------------------------------------------------*\
                       Class fvBoundaryMesh Declaration
\*---------------------------------------------------------------------------*/

class fvBoundaryMesh
:
    public fvPatchList
{
    // Private Data

        //- Reference to mesh
        const fvMesh& mesh_;


    // Private Member Functions

        //- Add fvPatches corresponding to the given polyBoundaryMesh
        void addPatches(const polyBoundaryMesh&);


        //- No copy construct
        fvBoundaryMesh(const fvBoundaryMesh&) = delete;

        //- No copy assignment
        void operator=(const fvBoundaryMesh&) = delete;


protected:

        //- Update boundary based on new polyBoundaryMesh
        void readUpdate(const polyBoundaryMesh&);


public:

    friend class fvMesh;


    // Constructors

        //- Construct zero size with mesh reference
        explicit fvBoundaryMesh(const fvMesh&);

        //- Construct from polyBoundaryMesh
        fvBoundaryMesh(const fvMesh&, const polyBoundaryMesh&);


    // Member Functions

        //- Return the mesh reference
        const fvMesh& mesh() const noexcept
        {
            return mesh_;
        }

        //- Return a list of faceCells for each patch
        UPtrList<const labelUList> faceCells() const;

        //- Return a list of pointers for each patch
        //- with only those pointing to interfaces being set
        lduInterfacePtrsList interfaces() const;

        //- Return (sorted) patch indices for all matches.
        //  A no-op (returns empty list) for an empty matcher
        labelList indices(const wordRe& matcher, const bool useGroups) const;

        //- Return (sorted) patch indices for all matches.
        //  A no-op (returns empty list) for an empty matcher
        labelList indices(const wordRes& matcher, const bool useGroups) const;

        //- Find patch index given a name
        //  A no-op (returns -1) for an empty patchName
        label findPatchID(const word& patchName) const;

        //- Correct patches after moving points
        void movePoints();


    // Member Operators

        //- Return const and non-const reference to fvPatch by index.
        using fvPatchList::operator[];

        //- Return const reference to fvPatch by name.
        const fvPatch& operator[](const word&) const;

        //- Return reference to fvPatch by name.
        fvPatch& operator[](const word&);


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
