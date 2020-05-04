/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::DictionaryBase

Description
    Base dictionary class templated on both the form of doubly-linked list
    it uses as well as the type it holds.

    The double templating allows for the instantiation of forms with or
    without storage management.

Note
    The IDLListType parameter should itself be a template but this confused
    gcc 2.95.2 so it has to be instantiated for T when an instantiation of
    DictionaryBase is requested

See also
    Dictionary and UDictionary

SourceFiles
    DictionaryBase.C
    DictionaryBaseIO.C

\*---------------------------------------------------------------------------*/

#ifndef DictionaryBase_H
#define DictionaryBase_H

#include "HashTable.H"
#include "wordList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations

template<class IDLListType, class T>
class DictionaryBase;

template<class IDLListType, class T>
Ostream& operator<<(Ostream&, const DictionaryBase<IDLListType, T>&);


/*---------------------------------------------------------------------------*\
                      Class DictionaryBase Declaration
\*---------------------------------------------------------------------------*/

template<class IDLListType, class T>
class DictionaryBase
:
    public IDLListType
{
protected:

    // Protected Data

        //- HashTable of the entries held on the IDLListType for quick lookup
        HashTable<T*> hashedTs_;


    // Protected Member Functions

        // Add the IDLListType entries into the HashTable
        void addEntries();


public:

    // Constructors

        //- Construct with given or default (128) table capacity
        explicit DictionaryBase(const label size = 128);

        //- Copy construct
        DictionaryBase(const DictionaryBase& dict);

        //- Construct from Istream using given Istream constructor class
        template<class INew>
        DictionaryBase(Istream& is, const INew& inew);

        //- Construct from Istream using default Istream constructor class
        DictionaryBase(Istream& is);


    // Member Functions

        // Search and lookup

            //- Search for given keyword
            bool found(const word& keyword) const;

            //- Find and return an entry, nullptr on failure.
            const T* cfind(const word& keyword) const;

            //- Find and return an entry, nullptr on failure.
            T* find(const word& keyword);

            //- Find and return entry, FatalError on failure.
            const T* lookup(const word& keyword) const;

            //- Find and return entry, FatalError on failure.
            T* lookup(const word& keyword);

            //- Return the table of contents (as a sorted list)
            wordList toc() const;

            //- Return the table of contents as a sorted list
            wordList sortedToc() const;

            //- Return table of contents sorted using the specified comparator
            template<class Compare>
            wordList sortedToc(const Compare& comp) const;


        // Editing

            //- Add at head of dictionary
            void insert(const word& keyword, T*);

            //- Add at tail of dictionary
            void append(const word& keyword, T*);

            //- Remove and return entry specified by keyword.
            //  Return nullptr if the keyword was not found.
            T* remove(const word& keyword);

            //- Clear the dictionary
            void clear();

            //- Transfer the contents of the argument into this DictionaryBase
            //  and annul the argument.
            void transfer(DictionaryBase<IDLListType, T>& dict);


    // Member Operators

        //- Copy assignment
        void operator=(const DictionaryBase&);

        //- Find and return entry
        const T* operator[](const word& key) const
        {
            return lookup(key);
        }

        //- Find and return entry
        T* operator[](const word& key)
        {
            return lookup(key);
        }


    // Ostream Operator

        friend Ostream& operator<< <IDLListType, T>
        (
            Ostream&,
            const DictionaryBase<IDLListType, T>&
        );


    // Housekeeping

        //- Deprecated(2020-03) use cfind()
        //
        //  \deprecated(2020-03) - use cfind() method
        FOAM_DEPRECATED_FOR(2020-03, "cfind() method")
        const T* lookupPtr(const word& keyword) const
        {
            return this->cfind(keyword);
        }

        //- Deprecated(2020-03) use find()
        //
        //  \deprecated(2020-03) - use find() method
        FOAM_DEPRECATED_FOR(2020-03, "find() method")
        T* lookupPtr(const word& keyword)
        {
            return this->find(keyword);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DictionaryBase.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
