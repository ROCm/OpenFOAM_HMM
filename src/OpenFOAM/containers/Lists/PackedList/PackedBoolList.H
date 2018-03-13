/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017-2018 OpenCFD Ltd.
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
    Foam::PackedBoolList

Description
    A bit-packed bool list.

    In addition to the obvious memory advantage over using a
    List\<bool\>, this class also provides a number of bit-like
    operations.

SourceFiles
    PackedBoolListI.H
    PackedBoolList.C

See also
    Foam::PackedList

\*---------------------------------------------------------------------------*/

#ifndef PackedBoolList_H
#define PackedBoolList_H

#include "PackedList.H"
#include "UIndirectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration
class PackedBoolList;

//- \typedef A List of PackedBoolList
typedef List<PackedBoolList> PackedBoolListList;

/*---------------------------------------------------------------------------*\
                       Class PackedBoolList Declaration
\*---------------------------------------------------------------------------*/

class PackedBoolList
:
    public PackedList<1>
{
    // Private Member Functions

        //- Preparation, resizing before a bitor operation
        //  returns true if the later result needs trimming
        bool bitorPrepare(const PackedBoolList& lst, label& maxPackLen);

        //- Set the listed indices. Return number of elements changed.
        //  Does auto-vivify for non-existent entries.
        template<class LabelListType>
        void setIndices(const LabelListType& indices);

        //- Unset the listed indices. Return number of elements changed.
        //  Never auto-vivify entries.
        template<class LabelListType>
        void unsetIndices(const LabelListType& indices);

public:

    // Constructors

        //- Construct null
        PackedBoolList() = default;

        //- Construct from Istream
        PackedBoolList(Istream& is);

        //- Construct with given size, initializes list to 0 (false)
        explicit inline PackedBoolList(const label size);

        //- Construct with given size and value for all elements
        inline PackedBoolList(const label size, const bool val);

        //- Copy construct
        inline PackedBoolList(const PackedBoolList& lst);

        //- Move construct
        inline PackedBoolList(PackedBoolList&& lst);

        //- Construct with given size and list of labels to set as true.
        inline PackedBoolList(const label size, const labelUList& indices);

        //- Construct with given size and list of labels to set as true.
        inline PackedBoolList(const label size, const labelUIndList& indices);

        //- Construct from a list of bools
        explicit inline PackedBoolList(const UList<bool>& lst);

        //- Construct from a list of labels
        //- using the labels as indices to indicate which bits are set
        explicit inline PackedBoolList(const labelUList& indices);

        //- Construct from a list of labels
        //- using the labels as indices to indicate which bits are set
        explicit inline PackedBoolList(const labelUIndList& indices);

        //- Clone
        inline autoPtr<PackedBoolList> clone() const;


    // Member Functions

    // Query

        //- Test value at specified position.
        //  \note Method name compatibility with std::bitset
        inline bool test(const label pos) const;


    // Access

        //- Single index/value assign
        using PackedList<1>::set;

        //- Single index unassign
        using PackedList<1>::unset;

        //- Set specified bits.
        void set(const PackedBoolList& lst);

        //- Unset specified bits.
        void unset(const PackedBoolList& lst);

        //- Return indices of the used (true) elements as a list of labels
        labelList used() const;


    // Edit

        //- Transfer the contents of the argument list into this list
        //- and annul the argument list.
        inline void transfer(PackedBoolList& lst);


    // Convenience Methods

        //- Set the listed indices. Return number of elements changed.
        //  Does auto-vivify for non-existent entries.
        void setMany(const labelUList& indices);

        //- Set the listed indices. Return number of elements changed.
        //  Does auto-vivify for non-existent entries.
        void setMany(const labelUIndList& indices);

        //- Unset the listed indices. Return number of elements changed.
        //  Never auto-vivify entries.
        void unsetMany(const labelUList& indices);

        //- Unset the listed indices. Return number of elements changed.
        //  Never auto-vivify entries.
        void unsetMany(const labelUIndList& indices);


    // Member Operators

        //- Assign all entries to the given value.
        inline void operator=(const bool val);

        //- Copy assignment
        inline void operator=(const PackedBoolList& lst);

        //- Move assignment
        inline void operator=(PackedBoolList&& lst);


    // Housekeeping

        //- No assignment from list. Use setMany for that.
        void operator=(const labelUList&) = delete;

        //- No assignment from list. Use setMany for that.
        void operator=(const labelUIndList&) = delete;


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PackedBoolListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
