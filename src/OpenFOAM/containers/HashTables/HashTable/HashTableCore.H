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
    Foam::HashTableCore

Description
    Template invariant parts of hash table implementation.

SourceFiles
    HashTableCoreI.H
    HashTableCore.C

\*---------------------------------------------------------------------------*/

#ifndef HashTableCore_H
#define HashTableCore_H

#include "label.H"
#include "uLabel.H"
#include "className.H"
#include "nullObject.H"
#include "zero.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class HashTableCore Declaration
\*---------------------------------------------------------------------------*/

//- Bits that are independent of HashTable template parameters.
struct HashTableCore
{
    //- Maximum allowable internal table size. Approximately labelMax/4
    static const label maxTableSize;

    //- Return a canonical (power-of-two) of the requested size.
    static label canonicalSize(const label requested_size);

    //- Construct null
    HashTableCore()
    {}

    //- Define template name and debug
    ClassName("HashTable");

    static_assert
    (
        sizeof(NullObject) >= sizeof(void*),
        "NullObject is too small to reinterpret_cast as HashTable::iterator"
    );


    //- Factory method to create a non-const iterator begin
    template<class IteratorType, class TableType>
    inline static IteratorType iterator_begin(TableType& table);

    //- Factory method to create a const iterator begin
    template<class IteratorType, class TableType>
    inline static IteratorType iterator_cbegin(const TableType& table);

    //- Factory method to return an iterator end
    //  Simply reinterprets a NullObject as a hash-table iterator.
    template<class IteratorType>
    inline static const IteratorType& iterator_end();

    //- Factory method to return an iterator cend
    //  Simply reinterprets a NullObject as a hash-table iterator.
    template<class IteratorType>
    inline static const IteratorType& iterator_cend();


    //- Factory class for creating a begin/end pair for any const iterator.
    template<class IteratorType, class TableType>
    class const_iterator_pair
    {
        label size_;
        IteratorType iter_;

    public:

        inline const_iterator_pair(const TableType& tbl);

        inline label size() const;
        inline bool empty() const;

        inline IteratorType begin() const;
        inline IteratorType cbegin() const;

        inline const IteratorType& end() const;
        inline const IteratorType& cend() const;
    };


protected:

    //- A static zero::null for dereferencing as a dummy HashSet element
    static zero::null zeroNullElement;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "HashTableCoreI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
