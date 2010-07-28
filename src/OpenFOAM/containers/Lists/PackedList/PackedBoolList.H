/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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
    Foam::PackedBoolList

Description
    A bit-packed bool list

SourceFiles
    PackedBoolListI.H
    PackedBoolList.C

SeeAlso
    Foam::PackedList

\*---------------------------------------------------------------------------*/

#ifndef PackedBoolList_H
#define PackedBoolList_H

#include "PackedList.H"
#include "UIndirectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class PackedBoolList Declaration
\*---------------------------------------------------------------------------*/

class PackedBoolList
:
    public PackedList<1>
{

    // Private Member Functions

        //- Modulo is everything that is not in lst
        //  This is equivalent to unsetting the bits specified in lst
        void modulo(const PackedList<1>&);

        //- Set the listed indices. Return number of elements changed.
        //  Does auto-vivify for non-existent entries.
        template<class LabelListType>
        label setFromIndices(const LabelListType& indices);

        //- Unset the listed indices. Return number of elements changed.
        //  Never auto-vivify entries.
        template<class LabelListType>
        label unsetFromIndices(const LabelListType& indices);


public:

    // Constructors

        //- Construct null
        inline PackedBoolList();

        //- Construct with given size, initializes list to 0
        explicit inline PackedBoolList(const label size);

        //- Construct with given size and value for all elements
        inline PackedBoolList(const label size, const bool val);

        //- Copy constructor
        inline PackedBoolList(const PackedBoolList&);

        //- Copy constructor
        explicit inline PackedBoolList(const PackedList<1>&);

        //- Construct by transferring the parameter contents
        inline PackedBoolList(const Xfer<PackedBoolList>&);

        //- Construct by transferring the parameter contents
        inline PackedBoolList(const Xfer<PackedList<1> >&);

        //- Construct from a list of bools
        explicit inline PackedBoolList(const UList<bool>&);

        //- Construct from a list of labels
        //  using the labels as indices to indicate which bits are set
        explicit inline PackedBoolList(const UList<label>& indices);

        //- Construct from a list of labels
        //  using the labels as indices to indicate which bits are set
        explicit inline PackedBoolList(const UIndirectList<label>& indices);

        //- Clone
        inline autoPtr<PackedBoolList> clone() const;


    // Member Functions

        // Access

            //- Return indices of the used (true) elements as a list of labels
            Xfer<labelList> used() const;

            using PackedList<1>::set;
            using PackedList<1>::unset;

            //- Set the listed indices. Return number of elements changed.
            //  Does auto-vivify for non-existent entries.
            label set(const UList<label>& indices);

            //- Set the listed indices. Return number of elements changed.
            //  Does auto-vivify for non-existent entries.
            label set(const UIndirectList<label>& indices);

            //- Unset the listed indices. Return number of elements changed.
            //  Never auto-vivify entries.
            label unset(const UList<label>& indices);

            //- Unset the listed indices. Return number of elements changed.
            //  Never auto-vivify entries.
            label unset(const UIndirectList<label>& indices);


        // Edit

            //- Transfer the contents of the argument list into this list
            //  and annul the argument list.
            inline void transfer(PackedBoolList&);

            //- Transfer the contents of the argument list into this list
            //  and annul the argument list.
            inline void transfer(PackedList<1>&);

            //- Transfer contents to the Xfer container
            inline Xfer<PackedBoolList> xfer();


    // Member Operators

            //- Assignment of all entries to the given value.
            inline PackedBoolList& operator=(const bool val);

            //- Assignment operator.
            inline PackedBoolList& operator=(const PackedBoolList&);

            //- Assignment operator.
            inline PackedBoolList& operator=(const PackedList<1>&);

            //- Assignment operator.
            PackedBoolList& operator=(const UList<bool>&);

            //- Assignment operator,
            //  using the labels as indices to indicate which bits are set
            PackedBoolList& operator=(const UList<label>& indices);

            //- Assignment operator,
            //  using the labels as indices to indicate which bits are set
            PackedBoolList& operator=(const UIndirectList<label>& indices);

            //- Complement operator
            inline PackedBoolList operator~() const;

            //- Or operator (lists may be dissimilar sizes)
            PackedBoolList& operator|=(const PackedList<1>&);

            //- Or operator (lists may be dissimilar sizes),
            //  using the labels as indices to indicate which bits are set
            inline PackedBoolList& operator|=(const UList<label>& indices);

            //- Or operator (lists may be dissimilar sizes),
            //  using the labels as indices to indicate which bits are set
            inline PackedBoolList& operator|=(const UIndirectList<label>&);

            //- And operator (lists may be dissimilar sizes)
            //  The result is trimmed to the smallest intersecting size
            PackedBoolList& operator&=(const PackedList<1>&);

            //- And operator (lists may be dissimilar sizes)
            //  Only retain unique entries
            PackedBoolList& operator^=(const PackedList<1>&);


            //- Synonymous with the or operator
            inline PackedBoolList& operator+=(const PackedList<1>&);

            //- Synonymous with the or operator
            inline PackedBoolList& operator+=(const UList<label>& indices);

            //- Synonymous with the or operator
            inline PackedBoolList& operator+=(const UIndirectList<label>&);

            //- Remove entries from this list
            inline PackedBoolList& operator-=(const PackedList<1>&);

            //- Remove entries from this list
            inline PackedBoolList& operator-=(const UList<label>& indices);

            //- Remove entries from this list
            inline PackedBoolList& operator-=(const UIndirectList<label>&);

};


// Global Operators

//- Combine lists
PackedBoolList operator|
(
    const PackedBoolList& lst1,
    const PackedBoolList& lst2
);


//- Intersect lists - the result is trimmed to the smallest intersecting size
PackedBoolList operator&
(
    const PackedBoolList& lst1,
    const PackedBoolList& lst2
);


//- Combine to form a unique list (xor)
PackedBoolList operator^
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
