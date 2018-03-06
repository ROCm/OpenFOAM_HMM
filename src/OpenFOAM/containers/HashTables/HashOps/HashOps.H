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

Namspace
    Foam::HashSetOps

Description
    Various operations for HashSet.

Namspace
    Foam::HashTableOps

Description
    Various operations for HashTable.

SourceFiles
    HashOps.H

\*---------------------------------------------------------------------------*/

#ifndef HashOps_H
#define HashOps_H

#include "HashSet.H"
#include "List.H"

namespace Foam
{

// Forward declarations

class PackedBoolList;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                       Namespace HashSetOps Declaration
\*---------------------------------------------------------------------------*/

namespace HashSetOps
{

//- Combine HashSet operation. Equivalent to 'a += b'
template<class Key=word, class Hash=string::hash>
struct plusEqOp
{
    typedef HashSet<Key, Hash> value_type;

    void operator()(value_type& a, const value_type& b) const
    {
        a += b;
    }
};


//- Convert a packed list of bits to a labelHashSet of the indices used.
//
//  \param selection the list for which a 'true' entry corresponds
//       to an index for be added to the labelHashSet
//
//  \return a labelHashSet of the selected indices
labelHashSet used(const PackedBoolList& select);


//- Convert a list of bools to a labelHashSet of the indices used.
//
//  \param selection the list for which a 'true' entry corresponds
//       to an index for be added to the labelHashSet
//
//  \return a labelHashSet of the selected indices
labelHashSet used(const UList<bool>& select);


//- Convert labels to a packed list of bits, with '1' for each
//- non-negative value and '0' for all others.
//
//  \param labels the list of indices.
//
//  \return a packed bit list of the selected indices
//
//  \note The operation discards any negative values since these are
//     invalid positions in the output list.
PackedBoolList bitset(const labelHashSet& labels);


//- Convert labels to a list of bools, with 'true' for each
//- non-negative value and 'false' for all others.
//
//  \param labels the list of indices.
//
//  \return a bool List of the selected indices
//
//  \note The operation discards any negative values since these are
//     invalid positions in the output list.
List<bool> bools(const labelHashSet& labels);


} // End namespace HashSetOps


/*---------------------------------------------------------------------------*\
                      Namespace HashTableOps Declaration
\*---------------------------------------------------------------------------*/

namespace HashTableOps
{

//- Combine HashTable operation. Equivalent to 'a += b'
template<class T, class Key=word, class Hash=string::hash>
struct plusEqOp
{
    typedef HashTable<T, Key, Hash> value_type;

    void operator()(value_type& a, const value_type& b) const
    {
        a += b;
    }
};

} // End namespace HashTableOps


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
