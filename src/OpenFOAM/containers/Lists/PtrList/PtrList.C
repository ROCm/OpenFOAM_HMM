/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2018 OpenCFD Ltd.
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

#include "PtrList.H"
#include "SLPtrList.H"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class T>
Foam::PtrList<T>::PtrList(const PtrList<T>& list)
:
    UPtrList<T>(list.size())
{
    const label len = this->size();

    for (label i=0; i<len; ++i)
    {
        this->ptrs_[i] = (list[i]).clone().ptr();
    }
}


template<class T>
template<class CloneArg>
Foam::PtrList<T>::PtrList(const PtrList<T>& list, const CloneArg& cloneArg)
:
    UPtrList<T>(list.size())
{
    const label len = this->size();

    for (label i=0; i<len; ++i)
    {
        this->ptrs_[i] = (list[i]).clone(cloneArg).ptr();
    }
}


template<class T>
Foam::PtrList<T>::PtrList(PtrList<T>& list, bool reuse)
:
    UPtrList<T>(list, reuse)
{
    if (!reuse)
    {
        const label len = this->size();

        for (label i=0; i<len; ++i)
        {
            this->ptrs_[i] = (list[i]).clone().ptr();
        }
    }
}


template<class T>
Foam::PtrList<T>::PtrList(const SLPtrList<T>& list)
:
    UPtrList<T>(list.size())
{
    if (list.size())
    {
        label i = 0;
        for (auto iter = list.cbegin(); iter != list.cend(); ++iter)
        {
            this->ptrs_[i++] = (*iter).clone().ptr();
        }
    }
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

template<class T>
Foam::PtrList<T>::~PtrList()
{
    const label len = this->size();

    // Free old pointers
    for (label i=0; i<len; ++i)
    {
        if (this->ptrs_[i])
        {
            delete this->ptrs_[i];
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
void Foam::PtrList<T>::clear()
{
    const label len = this->size();

    for (label i=0; i<len; ++i)
    {
        if (this->ptrs_[i])
        {
            delete this->ptrs_[i];
        }
    }

    UPtrList<T>::clear();
}


template<class T>
void Foam::PtrList<T>::resize(const label newLen)
{
    if (newLen <= 0)
    {
        clear();
        return;
    }

    const label oldLen = this->size();

    if (newLen != oldLen)
    {
        // Truncation frees old pointers
        for (label i=newLen; i<oldLen; ++i)
        {
            if (this->ptrs_[i])
            {
                delete this->ptrs_[i];
            }
        }

        // Any new elements are initialized to nullptr.
        this->ptrs_.resize(newLen, reinterpret_cast<T*>(0));
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
void Foam::PtrList<T>::operator=(const PtrList<T>& list)
{
    if (this == &list)
    {
        FatalErrorInFunction
            << "attempted assignment to self for type " << typeid(T).name()
            << abort(FatalError);
    }

    const label oldLen = this->size();
    const label newLen = list.size();

    // Truncate (frees old pointers) or extend the length
    resize(newLen);

    if (newLen < oldLen)
    {
        // Copy values for existing entries
        for (label i=0; i<newLen; ++i)
        {
            (*this)[i] = list[i];
        }
    }
    else
    {
        // Copy values for existing entries
        for (label i=0; i<oldLen; ++i)
        {
            (*this)[i] = list[i];
        }

        // Clone pointers for new entries
        for (label i=oldLen; i<newLen; ++i)
        {
            this->ptrs_[i] = (list[i]).clone().ptr();
        }
    }
}


// ************************************************************************* //
