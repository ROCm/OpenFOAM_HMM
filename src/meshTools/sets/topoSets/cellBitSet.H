/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::cellBitSet

Description
    A special purpose topoSet with the cell labels stored as a bitSet.
    It does not correspond to a cellSet either (no associated IOobject).

SourceFiles
    cellBitSet.C

\*---------------------------------------------------------------------------*/

#ifndef cellBitSet_H
#define cellBitSet_H

#include "topoSet.H"
#include "bitSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class cellBitSet Declaration
\*---------------------------------------------------------------------------*/

class cellBitSet
:
    public topoSet
{
    // Private data

        bitSet selected_;

public:

    //- Runtime type information
    TypeName("cellBitSet");


    // Constructors

        //- Construct with nCells elements, all elements unset
        explicit cellBitSet(const polyMesh& mesh);

        //- Construct with nCells elements, using initial val
        cellBitSet(const polyMesh& mesh, const bool val);


    //- Destructor
    virtual ~cellBitSet() = default;


    // Member Functions

        //- Return the bitSet
        const bitSet& addressing() const
        {
            return selected_;
        }

        //- Access the bitSet
        bitSet& addressing()
        {
            return selected_;
        }

        //- Set values to false, leaving the size untouched
        void reset()
        {
            selected_.reset();
        }

        //- Has the given index?
        virtual bool found(const label id) const;

        //- Set an index
        virtual bool set(const label id);

        //- Unset an index
        virtual bool unset(const label id);

        //- Set multiple indices
        virtual void set(const labelUList& labels);

        //- Unset multiple indices
        virtual void unset(const labelUList& labels);

        //- Invert contents.
        //  Insert all members [0,maxLen) which were not in set.
        virtual void invert(const label maxLen);

        //- Subset contents. Only elements present in both sets remain.
        virtual void subset(const topoSet& set);

        //- Add elements present in set.
        virtual void addSet(const topoSet& set);

        //- Subtract elements present in set.
        virtual void subtractSet(const topoSet& set);

        //- Sync cellBitSet across coupled patches.
        virtual void sync(const polyMesh& mesh)
        {}

        //- Return max index+1.
        virtual label maxSize(const polyMesh& mesh) const;

        //- Update any stored data for new labels.
        virtual void updateMesh(const mapPolyMesh& morphMap)
        {}

        //- Update any stored data for mesh redistribution.
        virtual void distribute(const mapDistributePolyMesh& map)
        {}

        //- Write maxLen items with label and coordinates.
        virtual void writeDebug
        (
            Ostream& os,
            const primitiveMesh& mesh,
            const label maxLen
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
