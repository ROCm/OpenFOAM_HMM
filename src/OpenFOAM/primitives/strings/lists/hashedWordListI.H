/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

inline void Foam::hashedWordList::rehash(const bool unique)
{
    if (unique)
    {
        uniq();
    }
    else
    {
        rehash();
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::hashedWordList::hashedWordList()
:
    List<word>(),
    indices_()
{}


inline Foam::hashedWordList::hashedWordList(const hashedWordList& lst)
:
    List<word>(static_cast<const UList<word>&>(lst))
{
    rehash();
}


inline Foam::hashedWordList::hashedWordList
(
    const UList<word>& lst,
    const bool removeDuplicates
)
:
    List<word>(lst)
{
    rehash(removeDuplicates);
}


inline Foam::hashedWordList::hashedWordList
(
    const Xfer<List<word>>& lst,
    const bool removeDuplicates
)
:
    List<word>(lst)
{
    rehash(removeDuplicates);
}


inline Foam::hashedWordList::hashedWordList(std::initializer_list<word> lst)
:
    List<word>(lst)
{
    rehash();
}


template<class AnyType, class AnyHash>
inline Foam::hashedWordList::hashedWordList
(
    const HashTable<AnyType, word, AnyHash>& tbl
)
:
    List<word>(tbl.size())
{
    label count = 0;
    for (auto iter = tbl.cbegin(); iter != tbl.cend(); ++iter)
    {
        List<word>::operator[](count++) = iter.key();
    }

    this->sort();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::hashedWordList::clear()
{
    List<word>::clear();
    indices_.clear();
}


inline void Foam::hashedWordList::append
(
    const word& name,
    const bool avoidDuplicates
)
{
    // name is either unique or we don't care about duplicates
    if (indices_.insert(name, size()) || !avoidDuplicates)
    {
        List<word>::append(name);
    }
}


inline const Foam::HashTable<Foam::label>&
Foam::hashedWordList::lookup() const
{
    return indices_;
}


inline bool Foam::hashedWordList::found(const word& name) const
{
    return indices_.found(name);
}


inline bool Foam::hashedWordList::contains(const word& name) const
{
    return indices_.found(name);
}


inline void Foam::hashedWordList::sort()
{
    Foam::sort(*this);
    rehash();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline const Foam::word& Foam::hashedWordList::operator[]
(
    const label index
) const
{
    return List<word>::operator[](index);
}


inline Foam::label Foam::hashedWordList::operator[](const word& name) const
{
    auto iter = indices_.find(name);

    if (iter.found())
    {
        return iter.object();
    }

    return -1;  // Not found (or not hashed?)
}


inline bool Foam::hashedWordList::operator()(const word& name) const
{
    return indices_.found(name);
}


inline void Foam::hashedWordList::operator=(const UList<word>& lst)
{
    List<word>::operator=(lst);
    rehash();
}


inline void Foam::hashedWordList::operator=(std::initializer_list<word> lst)
{
    List<word>::operator=(lst);
    rehash();
}


inline void Foam::hashedWordList::operator=(const hashedWordList& lst)
{
    operator=(static_cast<const UList<word>&>(lst));
}


// ************************************************************************* //
