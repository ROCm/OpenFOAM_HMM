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
    A wordList with hashed indices for additional fast lookup by name.

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

class hashedWordList;

// Forward declaration of friend functions and operators
Istream& operator>>(Istream& is, hashedWordList& lst);
Ostream& operator<<(Ostream& os, const hashedWordList& lst);


/*---------------------------------------------------------------------------*\
                       Class hashedWordList Declaration
\*---------------------------------------------------------------------------*/

class hashedWordList
:
    public List<word>
{
    // Private data

        //- Hash of words/indices
        mutable HashTable<label> indices_;


    // Private Member Functions

        //- Rebuild the lookup hash or make unique entries first.
        inline void rehash(const bool unique);

public:

    // Constructors

        //- Construct null
        inline hashedWordList();

        //- Copy constructor.
        inline hashedWordList(const hashedWordList& lst);

        //- Construct from list of words,
        //  optionally eliminating duplicates
        inline hashedWordList
        (
            const UList<word>& lst,
            const bool removeDuplicates=false
        );

        //- Construct by transferring the parameter contents,
        //  optionally eliminating duplicates
        inline hashedWordList
        (
            const Xfer<List<word>>& lst,
            const bool removeDuplicates=false
        );

        //- Construct from an initializer list
        inline hashedWordList(std::initializer_list<word> lst);

        //- Construct from the word keys of any HashTable, sorting immediately.
        //  This also handles a wordHashSet, which is derived from a HashTable.
        //  The result is similar to a HashTable::sortedToc.
        template<class AnyType, class AnyHash>
        explicit inline hashedWordList
        (
            const HashTable<AnyType, word, AnyHash>& tbl
        );

        //- Construct from number and list of words,
        //  optionally eliminating duplicates
        hashedWordList
        (
            const label count,
            const char** lst,
            const bool removeDuplicates=false
        );

        //- Construct from a nullptr-terminated list of words,
        //  optionally eliminating duplicates
        hashedWordList(const char** lst, const bool removeDuplicates=false);

        //- Construct from Istream
        hashedWordList(Istream& is);


    // Member Functions

        //- Clear the list, i.e. set size to zero.
        inline void clear();

        //- Append an element at the end of the list,
        //  optionally avoid append if it would be a duplicate entry
        inline void append(const word& name, const bool avoidDuplicates=false);

        //- Does the list contain the specified name
        inline bool found(const word& name) const;

        //- Does the list contain the specified name
        inline bool contains(const word& name) const;

        //- Return the hash of words/indices for inspection
        inline const HashTable<label>& lookup() const;

        //- Transfer the contents of the argument List into this list
        //  and annul the argument list,
        //  optionally eliminating duplicates
        void transfer(List<word>& lst, const bool removeDuplicates=false);

        //- Rebuild the lookup hash indices
        void rehash() const;

        //- Sort the list and rehash the indices
        inline void sort();

        //- Adjust the list if necessary to eliminate duplicate entries,
        //  and rehash the indices
        void uniq();


    // Member Operators

        //- Return name corresponding to specified index
        inline const word& operator[](const label index) const;

        //- Return index corresponding to specified name, or -1 on failure
        inline label operator[](const word& name) const;

        //- Does the list contain the specified name - same as found.
        //  Makes hashedWordList suitable as a unary predicate.
        inline bool operator()(const word& name) const;


      // Assignment

        //- Assignment operator from list of words. Rehashes the indices.
        inline void operator=(const UList<word>& lst);

        //- Assignment operator from initializer list. Rehashes the indices.
        inline void operator=(std::initializer_list<word> lst);

        //- Assignment operator. Rehashes the indices.
        inline void operator=(const hashedWordList& lst);


    // Istream operators

        friend Istream& operator>>(Istream& is, hashedWordList& lst);
        friend Ostream& operator<<(Ostream& os, const hashedWordList& lst);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "hashedWordListI.H"

#endif

// ************************************************************************* //
