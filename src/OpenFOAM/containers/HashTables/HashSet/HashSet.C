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

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Key, class Hash>
bool Foam::HashSet<Key, Hash>::operator==(const HashSet<Key, Hash>& rhs) const
{
    // Are all lhs elements in rhs?
    for (const_iterator iter = this->begin(); iter != this->end(); ++iter)
    {
        if (!rhs.found(iter.key()))
        {
            return false;
        }
    }

    // Are all rhs elements in lhs?
    for (const_iterator iter = rhs.begin(); iter != rhs.end(); ++iter)
    {
        if (!found(iter.key()))
        {
            return false;
        }
    }

    return true;
}


template<class Key, class Hash>
bool Foam::HashSet<Key, Hash>::operator!=(const HashSet<Key, Hash>& rhs) const
{
    return !(operator==(rhs));
}


template<class Key, class Hash>
void Foam::HashSet<Key, Hash>::operator+=(const HashSet<Key, Hash>& rhs)
{
    // Add in rhs elements into lhs
    for (const_iterator iter = rhs.begin(); iter != rhs.end(); ++iter)
    {
        insert(iter.key());
    }
}


template<class Key, class Hash>
void Foam::HashSet<Key, Hash>::operator-=(const HashSet<Key, Hash>& rhs)
{
    // Remove rhs elements from lhs
    for (const_iterator iter = rhs.begin(); iter != rhs.end(); ++iter)
    {
        erase(iter.key());
    }
}


template<class Key, class Hash>
void Foam::HashSet<Key, Hash>::operator&=(const HashSet<Key, Hash>& rhs)
{
    // Remove elements not found in rhs as well
    for (iterator iter = this->begin(); iter != this->end(); ++iter)
    {
        if (!rhs.found(iter.key()))
        {
            erase(iter);
        }
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
