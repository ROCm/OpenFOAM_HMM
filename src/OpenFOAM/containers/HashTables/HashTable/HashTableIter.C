/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, class Key, class Hash>
template<bool Const>
Foam::HashTable<T, Key, Hash>::Iterator<Const>::Iterator
(
    table_type* tbl,
    const Key& key
)
:
    entry_(nullptr),
    container_(tbl),
    index_(0)
{
    if (tbl->size())
    {
        const label index = container_->hashKeyIndex(key);

        for (node_type* ep = container_->table_[index]; ep; ep = ep->next_)
        {
            if (key == ep->key())
            {
                entry_ = ep;
                index_ = index;
                break;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

//
// Any changes here may need changes in the iterator increment() method
//
template<class T, class Key, class Hash>
bool Foam::HashTable<T, Key, Hash>::iterator_erase
(
    node_type*& entry,
    label& index
)
{
    // Safeguard against the following:
    // - empty table
    // - nullptr entry
    // - end iterator (which is also a nullptr)
    // - negative index from a previous erase. See comment below.
    if (!size_ || !entry || index < 0)
    {
        return false;
    }

    // Decrease count
    --size_;

    // The previous element in the singly linked list
    node_type* prev = nullptr;

    for (node_type* ep = table_[index]; ep; ep = ep->next_)
    {
        if (ep == entry)
        {
            break;
        }
        prev = ep;
    }

    if (prev)
    {
        // Had previous element in linked list - reposition to there
        prev->next_ = entry->next_;
        delete entry;
        entry = prev;

        return true;
    }

    // Was first element on linked list
    table_[index] = entry->next_;
    delete entry;

    // Assign any non-nullptr value so it doesn't look like end()
    entry = reinterpret_cast<node_type*>(this);

    // Mark the present index to continue and bring it back to the present
    // location with the next index.
    //
    // Save: (-index-1), which has no ambiguity for index 0.
    // Retrieve:  (-(index+1))

    index = (-index - 1);

    return true;
}


// ************************************************************************* //
