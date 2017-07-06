/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::dictionary

Description
    A list of keyword definitions, which are a keyword followed by a number
    of values (eg, words and numbers) or by a sub-dictionary.
    Since the dictionary format is used extensively throughout OpenFOAM for
    input/output files, there are many examples of its use.

    Dictionary keywords are a plain word or a pattern (regular expression).
    The general order for searching is as follows:
      - exact match
      - pattern match (in reverse order)
      - optional recursion into the enclosing (parent) dictionaries

    The dictionary class is the base class for IOdictionary and also serves
    as a bootstrap dictionary for the objectRegistry data dictionaries.

Note
    Within dictionaries, entries can be referenced by using the '$' syntax
    familiar from shell programming. A '.' separator is used when referencing
    sub-dictionary entries. Leading '.' prefixes can be used to specify
    an entry from a parent directory. A leading ':' or '^' prefix specifies
    starting from the top-level entry. For example,

    \verbatim
    key1        val1;
    key2        $key1;   // use key1 value from current scope
    key3        $.key1;  // use key1 value from current scope

    subdict1
    {
        key1        val1b;
        key2        $..key1; // use key1 value from parent
        subdict2
        {
            key2        val2;
            key3        $...key1; // use key1 value from grandparent
        }
    }

    key4        $^subdict1.subdict2.key3;  // lookup with absolute scoping
    \endverbatim

    It is also possible to use the '${}' syntax for clarity.

SourceFiles
    dictionary.C
    dictionaryIO.C

SeeAlso
    - Foam::entry
    - Foam::dictionaryEntry
    - Foam::primitiveEntry

\*---------------------------------------------------------------------------*/

#ifndef dictionary_H
#define dictionary_H

#include "entry.H"
#include "IDLList.H"
#include "DLList.H"
#include "fileName.H"
#include "ITstream.H"
#include "HashTable.H"
#include "wordList.H"
#include "className.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators
class regExp;
class dictionary;
class SHA1Digest;

Istream& operator>>(Istream& is, dictionary& dict);
Ostream& operator<<(Ostream& os, const dictionary& dict);

/*---------------------------------------------------------------------------*\
                        Class dictionaryName Declaration
\*---------------------------------------------------------------------------*/

//- Holds name for a dictionary
class dictionaryName
{
    // Private data

        fileName name_;


public:

    // Constructors

        //- Construct null
        dictionaryName()
        {}

        //- Construct as copy of the given fileName
        dictionaryName(const fileName& name)
        :
            name_(name)
        {}


    // Member functions

        //- Return the dictionary name
        const fileName& name() const
        {
            return name_;
        }

        //- Return the dictionary name
        fileName& name()
        {
            return name_;
        }

        //- Return the local dictionary name (final part of scoped name)
        const word dictName() const
        {
            const word scopedName = name_.name();
            const std::string::size_type i = scopedName.rfind('.');

            if (i == std::string::npos)
            {
                return scopedName;
            }
            else
            {
                return scopedName.substr(i+1);
            }
        }
};


/*---------------------------------------------------------------------------*\
                           Class dictionary Declaration
\*---------------------------------------------------------------------------*/

class dictionary
:
    public dictionaryName,
    public IDLList<entry>
{
    // Private data

        //- Report optional keywords and values if not present in dictionary
        static bool writeOptionalEntries;

        //- Parent dictionary
        const dictionary& parent_;

        //- HashTable of the entries held on the IDLList for quick lookup
        HashTable<entry*> hashedEntries_;

        //- Entries of matching patterns
        DLList<entry*> patterns_;

        //- Patterns as precompiled regular expressions
        DLList<autoPtr<regExp>> regexps_;


    // Typedefs

        //- The storage container
        typedef IDLList<entry> parent_type;

        typedef DLList<entry*>::iterator pattern_iterator;
        typedef DLList<entry*>::const_iterator pattern_const_iterator;
        typedef DLList<autoPtr<regExp>>::iterator regexp_iterator;
        typedef DLList<autoPtr<regExp>>::const_iterator regexp_const_iterator;


    // Private Member Functions

        //- Find and return an entry data stream pointer if present
        //  otherwise return nullptr. Allows scoping using '.'
        const entry* lookupScopedSubEntryPtr
        (
            const word& keyword,
            bool recursive,
            bool patternMatch
        ) const;


public:

    //- Declare friendship with the entry class for IO
    friend class entry;


    // Declare name of the class and its debug switch
    ClassName("dictionary");

    //- An empty dictionary, which is also the parent for all dictionaries
    static const dictionary null;


    // Constructors

        //- Construct top-level dictionary null
        dictionary();

        //- Construct top-level empty dictionary with given name
        explicit dictionary(const fileName& name);

        //- Construct given the entry name, parent dictionary and Istream,
        //  reading entries until lastEntry or EOF
        dictionary
        (
            const fileName& name,
            const dictionary& parentDict,
            Istream& is
        );

        //- Construct top-level dictionary from Istream,
        //  reading entries until EOF
        dictionary(Istream& is);

        //- Construct top-level dictionary from Istream,
        //  reading entries until EOF, optionally keeping the header
        dictionary(Istream& is, const bool keepHeader);

        //- Construct as copy given the parent dictionary
        dictionary(const dictionary& parentDict, const dictionary&);

        //- Construct top-level dictionary as copy
        dictionary(const dictionary& dict);

        //- Construct top-level dictionary as copy from pointer to dictionary.
        //  A null pointer is treated like an empty dictionary.
        dictionary(const dictionary* dictPtr);

        //- Construct by transferring parameter contents given parent dictionary
        dictionary(const dictionary& parentDict, const Xfer<dictionary>& dict);

        //- Construct top-level dictionary by transferring parameter contents
        dictionary(const Xfer<dictionary>& dict);

        //- Construct and return clone
        autoPtr<dictionary> clone() const;

        //- Construct top-level dictionary on freestore from Istream
        static autoPtr<dictionary> New(Istream& is);


    //- Destructor
    virtual ~dictionary();


    // Member functions

      // Access

        //- Return the parent dictionary
        const dictionary& parent() const
        {
            return parent_;
        }

        //- Return the top of the tree
        const dictionary& topDict() const;

        //- Return line number of first token in dictionary
        label startLineNumber() const;

        //- Return line number of last token in dictionary
        label endLineNumber() const;

        //- Return the SHA1 digest of the dictionary contents
        SHA1Digest digest() const;

        //- Return the dictionary as a list of tokens
        tokenList tokens() const;


        // Search and lookup

            //- Search dictionary for given keyword
            //  If recursive, search parent dictionaries
            //  If patternMatch, use regular expressions
            bool found
            (
                const word& keyword,
                bool recursive = false,
                bool patternMatch = true
            ) const;

            //- Find and return an entry data stream pointer if present
            //  otherwise return nullptr.
            //  If recursive, search parent dictionaries.
            //  If patternMatch, use regular expressions
            const entry* lookupEntryPtr
            (
                const word& keyword,
                bool recursive,
                bool patternMatch
            ) const;

            //- Find and return an entry data stream pointer for manipulation
            //  if present otherwise return nullptr.
            //  If recursive, search parent dictionaries.
            //  If patternMatch, use regular expressions.
            entry* lookupEntryPtr
            (
                const word& keyword,
                bool recursive,
                bool patternMatch
            );

            //- Find and return an entry data stream if present otherwise error.
            //  If recursive, search parent dictionaries.
            //  If patternMatch, use regular expressions.
            const entry& lookupEntry
            (
                const word& keyword,
                bool recursive,
                bool patternMatch
            ) const;

            //- Find and return an entry data stream
            //  If recursive, search parent dictionaries.
            //  If patternMatch, use regular expressions.
            ITstream& lookup
            (
                const word& keyword,
                bool recursive = false,
                bool patternMatch = true
            ) const;

            //- Find and return a T,
            //  if not found throw a fatal error.
            //  If recursive, search parent dictionaries.
            //  If patternMatch, use regular expressions.
            template<class T>
            T lookupType
            (
                const word&,
                bool recursive=false,
                bool patternMatch=true
            ) const;

            //- Find and return a T, or return the given default value
            //  If recursive, search parent dictionaries.
            //  If patternMatch, use regular expressions.
            template<class T>
            T lookupOrDefault
            (
                const word& keyword,
                const T& deflt,
                bool recursive = false,
                bool patternMatch = true
            ) const;

            //- Find and return a T, if not found return the given
            //  default value, and add to dictionary.
            //  If recursive, search parent dictionaries.
            //  If patternMatch, use regular expressions.
            template<class T>
            T lookupOrAddDefault
            (
                const word& keyword,
                const T& deflt,
                bool recursive = false,
                bool patternMatch = true
            );

            //- Find an entry if present, and assign to T
            //  Returns true if the entry was found.
            //  If recursive, search parent dictionaries.
            //  If patternMatch, use regular expressions.
            template<class T>
            bool readIfPresent
            (
                const word& keyword,
                T& val,
                bool recursive = false,
                bool patternMatch = true
            ) const;

            //- Find and return an entry data stream pointer if present
            //  otherwise return nullptr. Allows scoping using '.'.
            //  Special handling for ':' at start of keyword and '..'.
            const entry* lookupScopedEntryPtr
            (
                const word& keyword,
                bool recursive,
                bool patternMatch
            ) const;

            //- Check if entry is a sub-dictionary
            bool isDict(const word& keyword) const;

            //- Find and return a sub-dictionary pointer if present
            //  otherwise return nullptr.
            const dictionary* subDictPtr(const word& keyword) const;

            //- Find and return a sub-dictionary pointer if present
            //  otherwise return nullptr.
            dictionary* subDictPtr(const word& keyword);

            //- Find and return a sub-dictionary
            const dictionary& subDict(const word& keyword) const;

            //- Find and return a sub-dictionary for manipulation
            dictionary& subDict(const word& keyword);

            //- Find and return a sub-dictionary as a copy, or
            //  return an empty dictionary if the sub-dictionary does not exist
            dictionary subOrEmptyDict
            (
                const word& keyword,
                const bool mustRead = false
            ) const;

            //- Find and return a sub-dictionary if found
            //  otherwise return this dictionary
            const dictionary& optionalSubDict(const word&) const;

            //- Return the table of contents
            wordList toc() const;

            //- Return the sorted table of contents
            wordList sortedToc() const;

            //- Return the list of available keys or patterns
            List<keyType> keys(bool patterns = false) const;


        // Editing

            //- Substitute the given keyword prepended by '$' with the
            //  corresponding sub-dictionary entries
            bool substituteKeyword(const word& keyword, bool mergeEntry=false);

            //- Substitute the given scoped keyword prepended by '$' with the
            //  corresponding sub-dictionary entries
            bool substituteScopedKeyword
            (
                const word& keyword,
                bool mergeEntry=false
            );

            //- Add a new entry
            //  With the merge option, dictionaries are interwoven and
            //  primitive entries are overwritten
            bool add(entry* entryPtr, bool mergeEntry=false);

            //- Add an entry
            //  With the merge option, dictionaries are interwoven and
            //  primitive entries are overwritten
            void add(const entry& e, bool mergeEntry=false);

            //- Add a word entry
            //  optionally overwrite an existing entry
            void add(const keyType& k, const word& w, bool overwrite=false);

            //- Add a string entry
            //  optionally overwrite an existing entry
            void add(const keyType& k, const string& s, bool overwrite=false);

            //- Add a label entry
            //  optionally overwrite an existing entry
            void add(const keyType&, const label l, bool overwrite=false);

            //- Add a scalar entry
            //  optionally overwrite an existing entry
            void add(const keyType&, const scalar s, bool overwrite=false);

            //- Add a dictionary entry
            //  optionally merge with an existing sub-dictionary
            void add
            (
                const keyType& k,
                const dictionary& d,
                bool mergeEntry = false
            );

            //- Add a T entry
            //  optionally overwrite an existing entry
            template<class T>
            void add(const keyType& k, const T& t, bool overwrite=false);

            //- Assign a new entry, overwrite any existing entry
            void set(entry* entryPtr);

            //- Assign a new entry, overwrite any existing entry
            void set(const entry& e);

            //- Assign a dictionary entry, overwrite any existing entry
            void set(const keyType& k, const dictionary& d);

            //- Assign a T entry, overwrite any existing entry
            template<class T>
            void set(const keyType& k, const T& t);

            //- Remove an entry specified by keyword
            bool remove(const word& Keyword);

            //- Change the keyword for an entry,
            //  optionally forcing overwrite of an existing entry
            bool changeKeyword
            (
                const keyType& oldKeyword,
                const keyType& newKeyword,
                bool forceOverwrite=false
            );

            //- Merge entries from the given dictionary.
            //  Also merge sub-dictionaries as required.
            bool merge(const dictionary& dict);

            //- Clear the dictionary
            void clear();

            //- Transfer the contents of the argument and annul the argument.
            void transfer(dictionary& dict);

            //- Transfer contents to the Xfer container
            Xfer<dictionary> xfer();


        // Read

            //- Read dictionary from Istream
            bool read(Istream& is);

            //- Read dictionary from Istream, optionally keeping the header
            bool read(Istream& is, const bool keepHeader);


        // Write

            //- Write sub-dictionary with the keyword as its header
            void writeEntry(const keyType& keyword, Ostream& os) const;

            //- Write dictionary entries.
            //  Optionally with extra new line between entries for
            //  "top-level" dictionaries
            void writeEntries(Ostream& os, const bool extraNewLine=false) const;

            //- Write dictionary, normally with sub-dictionary formatting
            void write(Ostream& os, const bool subDict=true) const;


    // Member Operators

        //- Find and return entry
        ITstream& operator[](const word& keyword) const;

        void operator=(const dictionary& rhs);

        //- Include entries from the given dictionary.
        //  Warn, but do not overwrite existing entries.
        void operator+=(const dictionary& rhs);

        //- Conditionally include entries from the given dictionary.
        //  Do not overwrite existing entries.
        void operator|=(const dictionary& rhs);

        //- Unconditionally include entries from the given dictionary.
        //  Overwrite existing entries.
        void operator<<=(const dictionary& rhs);


    // IOstream operators

        //- Read dictionary from Istream
        friend Istream& operator>>(Istream& is, dictionary& dict);

        //- Write dictionary to Ostream
        friend Ostream& operator<<(Ostream& os, const dictionary& dict);
};


// Global Operators

//- Combine dictionaries.
//  Starting from the entries in dict1 and then including those from dict2.
//  Warn, but do not overwrite the entries from dict1.
dictionary operator+(const dictionary& dict1, const dictionary& dict2);

//- Combine dictionaries.
//  Starting from the entries in dict1 and then including those from dict2.
//  Do not overwrite the entries from dict1.
dictionary operator|(const dictionary& dict1, const dictionary& dict2);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "dictionaryTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
