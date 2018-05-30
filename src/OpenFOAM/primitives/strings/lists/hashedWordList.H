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

Class
    Foam::hashedWordList

Description
    A wordList with hashed named lookup, which can be faster in some
    situations than using the normal list find/found methods.

SourceFiles
    hashedWordListI.H
    hashedWordList.C

\*---------------------------------------------------------------------------*/

#ifndef hashedWordList_H
#define hashedWordList_H

#include "wordList.H"
#include "HashTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class hashedWordList;
inline Istream& operator>>(Istream& is, hashedWordList& lst);


/*---------------------------------------------------------------------------*\
                       Class hashedWordList Declaration
\*---------------------------------------------------------------------------*/

class hashedWordList
:
    public wordList
{
    // Private Data

        //- Lookup HashTable of words vs list-indices
        mutable HashTable<label> lookup_;

public:

    // Constructors

        //- Construct an empty list
        inline hashedWordList() = default;

        //- Copy construct.
        inline hashedWordList(const hashedWordList& lst);

        //- Move construct.
        inline hashedWordList(hashedWordList&& lst);

        //- Copy construct from list of words
        inline explicit hashedWordList(const wordUList& lst);

        //- Copy construct from list of words, eliminating duplicates
        inline hashedWordList(const wordUList& lst, bool unique);

        //- Move construct from list of words, optionally eliminating duplicates
        inline hashedWordList(wordList&& lst, bool unique=false);

        //- Construct from an initializer list
        inline hashedWordList(std::initializer_list<word> lst);

        //- Construct from the word keys of any HashTable, sorting immediately.
        //  This also handles a wordHashSet, which is derived from a HashTable.
        //  The result is similar to a HashTable::sortedToc.
        template<class AnyType, class AnyHash>
        inline explicit hashedWordList
        (
            const HashTable<AnyType, word, AnyHash>& tbl
        );

        //- Construct from number and list of words,
        //  optionally eliminating duplicates
        hashedWordList(const label len, const char** array, bool unique=false);

        //- Construct from a nullptr-terminated list of words,
        //  optionally eliminating duplicates
        hashedWordList(const char** array, bool unique=false);

        //- Construct from Istream
        inline hashedWordList(Istream& is);


    // Member Functions

        //- Clear the list, i.e. set size to zero.
        inline void clear();

        //- Append an element at the end of the list,
        //  optionally avoid append if it would be a duplicate entry
        inline void append(const word& name, bool unique=false);

        //- Does the list contain the specified name
        inline bool found(const word& name) const;

        //- Does the list contain the specified name
        inline bool contains(const word& name) const;

        //- Return the hash of words/indices for inspection
        inline const HashTable<label>& lookup() const;

        //- Swap contents
        inline void swap(hashedWordList& lst);

        //- Transfer contents of the argument into this list
        //- and annul the argument list, optionally eliminating duplicates
        inline void transfer(hashedWordList& lst);

        //- Transfer the contents of the argument List into this list
        //- and annul the argument list, optionally eliminating duplicates
        inline void transfer(wordList& lst, bool unique=false);

        //- Rebuild the lookup hash indices
        void rehash() const;

        //- Rebuild the lookup hash indices, or make unique entries first.
        inline void rehash(bool unique);

        //- Adjust the list (if needed) to eliminate duplicate entries,
        //- and rehash the indices
        void uniq();

        //- Sort the list and rehash the indices
        inline void sort();


    // Member Operators

        //- Return name corresponding to specified index.
        //  Fatal for out of range values.
        inline const word& operator[](const label index) const;

        //- Return index corresponding to specified name, or -1 on failure
        inline label operator[](const word& name) const;

        //- Check hashed values for the specified name - same as found.
        //  Can be used as a unary predicate.
        inline bool operator()(const word& name) const;


      // Assignment

        //- Copy assignment. Rehashes the indices.
        inline void operator=(const hashedWordList& lst);

        //- Copy assignment from list of words. Rehashes the indices.
        inline void operator=(const wordUList& lst);

        //- Copy assignment from initializer list. Rehashes the indices.
        inline void operator=(std::initializer_list<word> lst);

        //- Move assignment operator.
        inline void operator=(hashedWordList&& lst);

        //- Move assignment from list of words. Rehashes the indices.
        inline void operator=(wordList&& lst);

        //- Read from an input stream. Rehashes the indices.
        inline friend Istream& operator>>(Istream& is, hashedWordList& lst);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "hashedWordListI.H"

#endif

// ************************************************************************* //
