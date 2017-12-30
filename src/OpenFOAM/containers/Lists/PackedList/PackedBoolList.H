/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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
        bool bitorPrepare(const PackedList<1>& lst, label& maxPackLen);

        //- Set the listed indices. Return number of elements changed.
        //  Does auto-vivify for non-existent entries.
        template<class LabelListType>
        label setIndices(const LabelListType& indices);

        //- Unset the listed indices. Return number of elements changed.
        //  Never auto-vivify entries.
        template<class LabelListType>
        label unsetIndices(const LabelListType& indices);

        //- Subset with the listed indices. Return number of elements subsetted.
        template<class LabelListType>
        label subsetIndices(const LabelListType& indices);


public:

    // Constructors

        //- Construct null
        inline PackedBoolList();

        //- Construct from Istream
        PackedBoolList(Istream& is);

        //- Construct with given size, initializes list to 0 (false)
        explicit inline PackedBoolList(const label size);

        //- Construct with given size and value for all elements
        inline PackedBoolList(const label size, const bool val);

        //- Copy constructor
        inline PackedBoolList(const PackedBoolList& lst);

        //- Copy constructor
        explicit inline PackedBoolList(const PackedList<1>& lst);

        //- Construct by transferring the parameter contents
        inline PackedBoolList(const Xfer<PackedBoolList>& lst);

        //- Construct by transferring the parameter contents
        inline PackedBoolList(const Xfer<PackedList<1>>& lst);

        //- Construct with given size and list of labels to set as true.
        inline PackedBoolList(const label size, const labelUList& indices);

        //- Construct with given size and list of labels to set as true.
        inline PackedBoolList(const label size, const labelUIndList& indices);

        //- Construct from a list of bools
        explicit inline PackedBoolList(const UList<bool>& lst);

        //- Construct from a list of labels
        //  using the labels as indices to indicate which bits are set
        explicit inline PackedBoolList(const labelUList& indices);

        //- Construct from a list of labels
        //  using the labels as indices to indicate which bits are set
        explicit inline PackedBoolList(const labelUIndList& indices);

        //- Clone
        inline autoPtr<PackedBoolList> clone() const;


    // Member Functions

      // Access

        using PackedList<1>::set;
        using PackedList<1>::unset;

        //- Set specified bits.
        void set(const PackedList<1>& lst);

        //- Set the listed indices. Return number of elements changed.
        //  Does auto-vivify for non-existent entries.
        label set(const labelUList& indices);

        //- Set the listed indices. Return number of elements changed.
        //  Does auto-vivify for non-existent entries.
        label set(const labelUIndList& indices);

        //- Unset specified bits.
        void unset(const PackedList<1>& lst);

        //- Unset the listed indices. Return number of elements changed.
        //  Never auto-vivify entries.
        label unset(const labelUList& indices);

        //- Unset the listed indices. Return number of elements changed.
        //  Never auto-vivify entries.
        label unset(const labelUIndList& indices);

        //- Subset with the specified list.
        void subset(const PackedList<1>& lst);

        //- Subset with the listed indices.
        //  Return number of elements subsetted.
        label subset(const labelUList& indices);

        //- Subset with the listed indices.
        //  Return number of elements subsetted.
        label subset(const labelUIndList& indices);


        //- Return indices of the used (true) elements as a list of labels
        Xfer<labelList> used() const;


      // Edit

        //- Transfer the contents of the argument list into this list
        //  and annul the argument list.
        inline void transfer(PackedBoolList& lst);

        //- Transfer the contents of the argument list into this list
        //  and annul the argument list.
        inline void transfer(PackedList<1>& lst);

        //- Transfer contents to the Xfer container
        inline Xfer<PackedBoolList> xfer();


    // Member Operators

        //- Assignment of all entries to the given value.
        inline void operator=(const bool val);

        //- Assignment operator.
        inline void operator=(const PackedBoolList& lst);

        //- Assignment operator.
        inline void operator=(const PackedList<1>& lst);

        //- Assignment operator.
        void operator=(const UList<bool>& lst);

        //- Assignment operator,
        //  using the labels as indices to indicate which bits are set
        inline void operator=(const labelUList& indices);

        //- Assignment operator,
        //  using the labels as indices to indicate which bits are set
        inline void operator=(const labelUIndList& indices);

        //- Complement operator
        inline PackedBoolList operator~() const;

        //- And operator (lists may be dissimilar sizes)
        inline PackedBoolList& operator&=(const PackedList<1>& lst);

        //- And operator (lists may be dissimilar sizes)
        //  using the labels as indices to indicate which bits are set
        inline PackedBoolList& operator&=(const labelUList& indices);

        //- And operator (lists may be dissimilar sizes)
        //  using the labels as indices to indicate which bits are set
        inline PackedBoolList& operator&=(const labelUIndList& indices);

        //- Xor operator (lists may be dissimilar sizes)
        //  Retains unique entries
        PackedBoolList& operator^=(const PackedList<1>& lst);

        //- Or operator (lists may be dissimilar sizes)
        inline PackedBoolList& operator|=(const PackedList<1>& lst);

        //- Or operator (lists may be dissimilar sizes),
        //  using the labels as indices to indicate which bits are set
        inline PackedBoolList& operator|=(const labelUList& indices);

        //- Or operator (lists may be dissimilar sizes),
        //  using the labels as indices to indicate which bits are set
        inline PackedBoolList& operator|=(const labelUIndList& indices);


        //- Add entries to this list, synonymous with the or operator
        inline PackedBoolList& operator+=(const PackedList<1>& lst);

        //- Add entries to this list, synonymous with the or operator
        inline PackedBoolList& operator+=(const labelUList& indices);

        //- Add entries to this list, synonymous with the or operator
        inline PackedBoolList& operator+=(const labelUIndList& indices);

        //- Remove entries from this list - unset the specified bits
        inline PackedBoolList& operator-=(const PackedList<1>& lst);

        //- Remove entries from this list - unset the specified bits
        inline PackedBoolList& operator-=(const labelUList& indices);

        //- Remove entries from this list - unset the specified bits
        inline PackedBoolList& operator-=(const labelUIndList& indices);
};


// Global Operators

//- Intersect lists - the result is trimmed to the smallest intersecting size
PackedBoolList operator&
(
    const PackedBoolList& lst1,
    const PackedBoolList& lst2
);


//- Combine to form a unique list (xor)
//  The result is trimmed to the smallest intersecting size
PackedBoolList operator^
(
    const PackedBoolList& lst1,
    const PackedBoolList& lst2
);


//- Combine lists
PackedBoolList operator|
(
    const PackedBoolList& lst1,
    const PackedBoolList& lst2
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "PackedBoolListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
