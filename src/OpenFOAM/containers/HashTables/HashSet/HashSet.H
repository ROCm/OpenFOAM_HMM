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

Class
    Foam::HashSet

Description
    A HashTable with keys but without contents.

Typedef
    Foam::wordHashSet

Description
    A HashSet with (the default) word keys.

Typedef
    Foam::labelHashSet

Description
    A HashSet with label keys.

\*---------------------------------------------------------------------------*/

#ifndef HashSet_H
#define HashSet_H

#include "HashTable.H"
#include "empty.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class HashSet Declaration
\*---------------------------------------------------------------------------*/

template<class Key=word, class Hash=string::hash>
class HashSet
:
    public HashTable<empty, Key, Hash>
{

public:

    typedef typename HashTable<empty, Key, Hash>::iterator iterator;
    typedef typename HashTable<empty, Key, Hash>::const_iterator const_iterator;


    // Constructors

        //- Construct given initial size
        HashSet(label size = 100)
        :
            HashTable<empty, Key, Hash>(size)
        {}

        //- Construct from Istream
        HashSet(Istream& is)
        :
            HashTable<empty, Key, Hash>(is)
        {}

        //- Construct from UList of Key
        HashSet(const UList<Key>& lst)
        :
            HashTable<empty, Key, Hash>(2*lst.size())
        {
            forAll(lst, i)
            {
                insert(lst[i]);
            }
        }

        //- Construct as copy
        HashSet(const HashSet<Key, Hash>& hs)
        :
            HashTable<empty, Key, Hash>(hs)
        {}

        //- Construct by transferring the parameter contents
        HashSet(const xfer<HashSet<Key, Hash> >& hs)
        :
            HashTable<empty, Key, Hash>(hs)
        {}


    // Member Functions

        // Edit

            //- Insert a new entry
            bool insert(const Key& key)
            {
                return HashTable<empty, Key, Hash>::insert(key, empty());
            }

            //- Same as insert (cannot overwrite empty content)
            bool set(const Key& key)
            {
                return HashTable<empty, Key, Hash>::insert(key, empty());
            }


    // Member Operators

        //- Equality. Two hashtables are equal when their contents are equal.
        //  Independent of table size or order.
        bool operator==(const HashSet<Key, Hash>&) const;

        //- The opposite of the equality operation.
        bool operator!=(const HashSet<Key, Hash>&) const;


        //- Add entries listed in the given HashSet to this HashSet
        void operator+=(const HashSet<Key, Hash>&);

        //- Remove entries listed in the given HashSet from this HashSet
        void operator-=(const HashSet<Key, Hash>&);

        //- Only retain entries found in both HashSets
        void operator&=(const HashSet<Key, Hash>&);

};


//- A HashSet with word keys.
typedef HashSet<> wordHashSet;

//- A HashSet with label keys.
typedef HashSet<label, Hash<label> > labelHashSet;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "HashSet.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
