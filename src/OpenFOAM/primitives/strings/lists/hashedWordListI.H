/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2018 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::hashedWordList::hashedWordList(const hashedWordList& lst)
:
    hashedWordList(static_cast<const wordUList&>(lst), false)
{}


inline Foam::hashedWordList::hashedWordList(hashedWordList&& lst)
:
    wordList(std::move(static_cast<wordList&>(lst))),
    lookup_(std::move(lst.lookup_))
{}


inline Foam::hashedWordList::hashedWordList(const wordUList& lst)
:
    hashedWordList(lst, false)
{}


inline Foam::hashedWordList::hashedWordList(const wordUList& lst, bool unique)
:
    wordList(lst)
{
    rehash(unique);
}


inline Foam::hashedWordList::hashedWordList(wordList&& lst, bool unique)
:
    wordList(std::move(lst))
{
    rehash(unique);
}


inline Foam::hashedWordList::hashedWordList(std::initializer_list<word> lst)
:
    wordList(lst)
{
    rehash();
}


template<class AnyType, class AnyHash>
inline Foam::hashedWordList::hashedWordList
(
    const HashTable<AnyType, word, AnyHash>& tbl
)
:
    wordList(tbl.size())
{
    wordList& list = *this;

    label count = 0;
    for (auto iter = tbl.cbegin(); iter != tbl.cend(); ++iter)
    {
        list[count++] = iter.key();
    }

    this->sort();
}


inline Foam::hashedWordList::hashedWordList(Istream& is)
{
    is  >> *this;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::hashedWordList::clear()
{
    wordList::clear();
    lookup_.clear();
}


inline void Foam::hashedWordList::append(const word& name, bool unique)
{
    // name is either unique or we don't care about duplicates
    if (lookup_.insert(name, size()) || !unique)
    {
        wordList::append(name);
    }
}


inline const Foam::HashTable<Foam::label>&
Foam::hashedWordList::lookup() const
{
    const label lenList = wordList::size();
    const label lenHash = lookup_.size();

    if ((lenList < lenHash) || (lenList && !lenHash))
    {
        rehash();   // Was somehow out of sync
    }

    return lookup_;
}


inline bool Foam::hashedWordList::found(const word& name) const
{
    return lookup().found(name);
}


inline bool Foam::hashedWordList::contains(const word& name) const
{
    return lookup().found(name);
}


inline void Foam::hashedWordList::swap(hashedWordList& lst)
{
    wordList::swap(static_cast<wordList&>(lst));
    lookup_.swap(lst.lookup_);
}


inline void Foam::hashedWordList::transfer
(
    hashedWordList& lst
)
{
    wordList::transfer(static_cast<wordList&>(lst));
    lookup_.transfer(lst.lookup_);
}


inline void Foam::hashedWordList::transfer(wordList& lst, bool unique)
{
    wordList::transfer(lst);
    rehash(unique);
}


inline void Foam::hashedWordList::rehash(bool unique)
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
    return wordList::operator[](index);
}


inline Foam::label Foam::hashedWordList::operator[](const word& name) const
{
    return lookup_.lookup(name, -1);  // -1 = not found or not hashed
}


inline bool Foam::hashedWordList::operator()(const word& name) const
{
    return lookup_.found(name);
}


inline void Foam::hashedWordList::operator=(const hashedWordList& lst)
{
    wordList::operator=(lst);
    rehash();
}


inline void Foam::hashedWordList::operator=(const wordUList& lst)
{
    wordList::operator=(lst);
    rehash();
}


inline void Foam::hashedWordList::operator=(std::initializer_list<word> lst)
{
    wordList::operator=(lst);
    rehash();
}


inline void Foam::hashedWordList::operator=(hashedWordList&& lst)
{
    wordList::transfer(static_cast<wordList&>(lst));
    lookup_ = std::move(lst.lookup_);
}


inline void Foam::hashedWordList::operator=(wordList&& lst)
{
    wordList::transfer(lst);
    rehash();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

inline Foam::Istream& Foam::operator>>(Istream& is, hashedWordList& lst)
{
    is  >> static_cast<wordList&>(lst);
    lst.rehash();

    return is;
}


// ************************************************************************* //
