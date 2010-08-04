/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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
    Foam::polyPatchID

Description
    A class holds the data needed to identify a patch in a dynamic mesh.

    The patch is identified by name and its index in the boundary mesh
    is updated if the mesh has changed.

\*---------------------------------------------------------------------------*/

#ifndef polyPatchID_H
#define polyPatchID_H

#include "keyType.H"
#include "labelList.H"
#include "polyBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class polyPatchID;
Ostream& operator<<(Ostream&, const polyPatchID&);


/*---------------------------------------------------------------------------*\
                         Class polyPatchID Declaration
\*---------------------------------------------------------------------------*/

class polyPatchID
{
    // Private data

        //- Patch name
        keyType key_;

        //- Patch indices
        labelList indices_;


public:

    // Constructors

        //- Construct from name
        polyPatchID(const keyType& key, const polyBoundaryMesh& bm)
        :
            key_(key),
            index_(bm.findIndices(key_))
        {}

        //- Construct from Istream
        polyPatchID(Istream& is, const polyBoundaryMesh& bm)
        :
            key_(is),
            index_(bm.findIndices(key_))
        {}


    // Destructor - default


    // Member Functions

        // Access

            //- Return name
            const keyType& name() const
            {
                return key_;
            }

            //- Return indices of matching patches
            const labelList& indices() const
            {
                return indices_;
            }

            //- Return index of first matching patch
            label index() const
            {
                return indices_.empty() ? -1 : indices_[0];
            }

            //- Has the patch been found
            bool active() const
            {
                return !indices_.empty();
            }


        // Edit

            //- Update
            void update(const polyBoundaryMesh& bm)
            {
                indices_ = bm.findIndices(key_);
            }


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const polyPatchID&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Ostream& operator<<(Ostream& os, const polyPatchID& p)
{
    os  << token::BEGIN_LIST
        << p.name() << token::SPACE << p.index()
        << token::END_LIST;

    // Check state of Ostream
    os.check("Ostream& operator<<(Ostream&, const polyPatchID&)");

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
