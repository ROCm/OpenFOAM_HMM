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

InNamspace
    Foam

Description
    Various functions to operate on HashTables.

SourceFiles
    HashTableOps.H

\*---------------------------------------------------------------------------*/

#ifndef HashTableOps_H
#define HashTableOps_H

#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Combine HashSet operation. Equivalent to 'a += b'
template<class Key=word, class Hash=string::hash>
struct HashSetPlusEqOp
{
    typedef HashSet<Key, Hash> value_type;

    void operator()(value_type& a, const value_type& b) const
    {
        a += b;
    }
};


//- Combine HashTable operation. Equivalent to 'a += b'
template<class T, class Key=word, class Hash=string::hash>
struct HashTablePlusEqOp
{
    typedef HashTable<T, Key, Hash> value_type;

    void operator()(value_type& a, const value_type& b) const
    {
        a += b;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
