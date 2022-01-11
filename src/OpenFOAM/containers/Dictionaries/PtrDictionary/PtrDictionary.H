/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

Class
    Foam::PtrDictionary

Description
    Template dictionary class which manages the storage associated with it.

    It is derived from DictionaryBase instantiated on a memory managed form of
    intrusive doubly-linked list of \<T\>.

\*---------------------------------------------------------------------------*/

#ifndef Foam_PtrDictionary_H
#define Foam_PtrDictionary_H

#include "DictionaryBase.H"
#include "DLPtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class PtrDictionary Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class PtrDictionary
:
    public DictionaryBase<DLPtrList<T>, T>
{
public:

    // Constructors

        //- Construct given initial table size
        explicit PtrDictionary(const label size = 128)
        :
            DictionaryBase<DLPtrList<T>, T>(size)
        {}

        //- Copy construct
        PtrDictionary(const PtrDictionary& dict)
        :
            DictionaryBase<DLPtrList<T>, T>(dict)
        {}

        //- Construct from Istream using given Istream constructor class
        template<class INew>
        PtrDictionary(Istream& is, const INew& inew)
        :
            DictionaryBase<DLPtrList<T>, T>(is, inew)
        {}

        //- Construct from Istream
        explicit PtrDictionary(Istream& is)
        :
            DictionaryBase<DLPtrList<T>, T>(is)
        {}


    // Member Operators

        //- Find and return entry
        const T& operator[](const word& key) const
        {
            return *DictionaryBase<DLPtrList<T>, T>::operator[](key);
        }

        //- Find and return entry
        T& operator[](const word& key)
        {
            return *DictionaryBase<DLPtrList<T>, T>::operator[](key);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
