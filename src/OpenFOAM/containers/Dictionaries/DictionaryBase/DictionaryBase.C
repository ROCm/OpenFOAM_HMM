/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "DictionaryBase.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class IDLListType, class T>
void Foam::DictionaryBase<IDLListType, T>::addEntries()
{
    for (auto iter = this->begin(); iter != this->end(); ++iter)
    {
        this->hashedTs_.insert((*iter).keyword(), &(*iter));
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class IDLListType, class T>
Foam::DictionaryBase<IDLListType, T>::DictionaryBase(const label size)
:
    hashedTs_(size)
{}


template<class IDLListType, class T>
Foam::DictionaryBase<IDLListType, T>::DictionaryBase
(
    const DictionaryBase& dict
)
:
    IDLListType(dict)
{
    addEntries();
}


template<class IDLListType, class T>
template<class INew>
Foam::DictionaryBase<IDLListType, T>::DictionaryBase
(
    Istream& is,
    const INew& iNew
)
:
    IDLListType(is, iNew)
{
    addEntries();
}


template<class IDLListType, class T>
Foam::DictionaryBase<IDLListType, T>::DictionaryBase(Istream& is)
:
    IDLListType(is)
{
    addEntries();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class IDLListType, class T>
bool Foam::DictionaryBase<IDLListType, T>::found(const word& keyword) const
{
    return hashedTs_.found(keyword);
}


template<class IDLListType, class T>
const T* Foam::DictionaryBase<IDLListType, T>::cfind
(
    const word& keyword
) const
{
    const auto iter = hashedTs_.cfind(keyword);

    if (iter.found())
    {
        return *iter;
    }

    return nullptr;
}


template<class IDLListType, class T>
T* Foam::DictionaryBase<IDLListType, T>::find(const word& keyword)
{
    auto iter = hashedTs_.find(keyword);

    if (iter.found())
    {
        return *iter;
    }

    return nullptr;
}


template<class IDLListType, class T>
const T* Foam::DictionaryBase<IDLListType, T>::lookup(const word& keyword) const
{
    const auto iter = hashedTs_.cfind(keyword);

    if (!iter.found())
    {
        FatalErrorInFunction
            << "'" << keyword << "' not found"
            << exit(FatalError);
    }

    return *iter;
}


template<class IDLListType, class T>
T* Foam::DictionaryBase<IDLListType, T>::lookup(const word& keyword)
{
    auto iter = hashedTs_.find(keyword);

    if (!iter.found())
    {
        FatalErrorInFunction
            << "'" << keyword << "' not found"
            << exit(FatalError);
    }

    return *iter;
}


template<class IDLListType, class T>
Foam::wordList Foam::DictionaryBase<IDLListType, T>::toc() const
{
    // Cannot rely on the items themselves having a keyword() method
    // so simply return the toc() from the hashed entries
    // Make it sorted, since anything else would have no meaning.
    return hashedTs_.sortedToc();
}


template<class IDLListType, class T>
Foam::wordList Foam::DictionaryBase<IDLListType, T>::sortedToc() const
{
    return hashedTs_.sortedToc();
}


template<class IDLListType, class T>
template<class Compare>
Foam::wordList Foam::DictionaryBase<IDLListType, T>::sortedToc
(
    const Compare& comp
) const
{
    return hashedTs_.sortedToc(comp);
}


template<class IDLListType, class T>
void Foam::DictionaryBase<IDLListType, T>::insert(const word& keyword, T* tPtr)
{
    // NOTE: we should probably check that HashTable::insert actually worked
    hashedTs_.insert(keyword, tPtr);
    IDLListType::insert(tPtr);
}


template<class IDLListType, class T>
void Foam::DictionaryBase<IDLListType, T>::append(const word& keyword, T* tPtr)
{
    // NOTE: we should probably check that HashTable::insert actually worked
    hashedTs_.insert(keyword, tPtr);
    IDLListType::append(tPtr);
}


template<class IDLListType, class T>
T* Foam::DictionaryBase<IDLListType, T>::remove(const word& keyword)
{
    auto iter = hashedTs_.find(keyword);

    if (iter.found())
    {
        T* ptr = IDLListType::remove(iter());
        hashedTs_.erase(iter);
        return ptr;
    }

    return nullptr;
}


template<class IDLListType, class T>
void Foam::DictionaryBase<IDLListType, T>::clear()
{
    IDLListType::clear();
    hashedTs_.clear();
}


template<class IDLListType, class T>
void Foam::DictionaryBase<IDLListType, T>::transfer
(
    DictionaryBase<IDLListType, T>& dict
)
{
    if (this == &dict)
    {
        return;  // Self-assignment is a no-op
    }

    IDLListType::transfer(dict);
    hashedTs_.transfer(dict.hashedTs_);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class IDLListType, class T>
void Foam::DictionaryBase<IDLListType, T>::operator=
(
    const DictionaryBase<IDLListType, T>& dict
)
{
    if (this == &dict)
    {
        return;  // Self-assignment is a no-op
    }

    IDLListType::operator=(dict);
    this->hashedTs_.clear();
    this->addEntries();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "DictionaryBaseIO.C"

// ************************************************************************* //
