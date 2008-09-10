/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#ifndef HashSet_C
#define HashSet_C

#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Key, class Hash>
bool HashSet<Key, Hash>::operator==(const HashSet<Key, Hash>& ht) const
{
    const HashTable<empty, Key, Hash>& a = *this;

    // Are all my elements in ht?
    for
    (
        typename HashTable<empty, Key, Hash>::const_iterator iter = a.begin();
        iter != a.end();
        ++iter
    )
    {
        if (!ht.found(iter.key()))
        {
            return false;
        }
    }

    // Are all ht elements in me?
    for
    (
        typename HashTable<empty, Key, Hash>::const_iterator iter = ht.begin();
        iter != ht.end();
        ++iter
    )
    {
        if (!found(iter.key()))
        {
            return false;
        }
    }
    return true;
}


template<class Key, class Hash>
bool HashSet<Key, Hash>::operator!=(const HashSet<Key, Hash>& ht) const
{
    return !(operator==(ht));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
