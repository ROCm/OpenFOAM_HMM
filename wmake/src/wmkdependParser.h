/*---------------------------------*- C++ -*---------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

@file wmkdependParser.atg

Description
    An attributed Coco/R grammar to parse C/C++, Fortran and Java files
    for include and import statements.

SourceFiles
    generated

\*---------------------------------------------------------------------------*/
// This file was generated with Coco/R C++ (7 Feb 2010)
// http://www.ssw.uni-linz.ac.at/coco/
// with these defines:
//     - FORCE_UTF8


#ifndef COCO_wmkdependPARSER_H__
#define COCO_wmkdependPARSER_H__

#include <iostream>
#include <string>
#include <list>

//! @brief A simple HashTable implementation
/**
 * @note This hash table is only vaguely STL-like. In accordance with
 * its present purpose, this hash table only supports a constIterator
 * and no deletions. For simplicity, the constIterator increment is
 * simply via a next() method. Instead of comparing to an end value,
 * the constIterator valid() method is used.
 * For example,
 * @code
 *    for
 *    (
 *         HashTable<foo>::constIterator iter = myHash.begin();
 *         iter.valid();
 *         iter.next()
 *    )
 *    {
 *        std::cerr<< "key: " << iter.key() << "\n";
 *    }
 * @endcode
 *
 */
class StringHashSet
{
    //! An entry within the HashTable
    struct hashedEntry
    {
        const std::string key_;   //<! The lookup key
        hashedEntry *next_;       //<! Pointer to next hashedEntry in sub-list

        hashedEntry(const std::string& key, hashedEntry *next=0)
        :
            key_(key), next_(next)
        {}
    };

    const int size_;   //<! fixed HashTable size
    hashedEntry** table_;

public:

    //! Construct with a default size
    StringHashSet(int size = 500)
    :
        size_(size),
        table_(new hashedEntry*[size_])
    {
        memset(table_, 0, size_ * sizeof(hashedEntry*));
    }

    //! Destructor
    ~StringHashSet()
    {
        for (int hashIdx = 0; hashIdx < size_; ++hashIdx)
        {
            hashedEntry* ep = table_[hashIdx];
            while (ep)
            {
                hashedEntry* del = ep;
                ep = ep->next_;
                delete del;
            }
        }
        delete[] table_;
        table_ = 0;
    }

    //! Return hash index for lookup name in hash table
    bool hashKeyIndex(const std::string& name) const
    {
        int hashIdx = 0;

        // calculate hash index
        for
        (
            std::string::const_iterator iter = name.begin();
            iter != name.end();
            ++iter
        )
        {
            hashIdx = hashIdx << 1 ^ *iter;
        }

        if (hashIdx < 0)
        {
            hashIdx = -hashIdx;
        }

        return hashIdx % size_;
    }


    //! Return true if name is found in hash table
    bool found(const std::string& name) const
    {
        const int hashIdx = hashKeyIndex(name);

        for (hashedEntry* ep = table_[hashIdx]; ep; ep = ep->next_)
        {
            if (name == ep->key_)
            {
                // found
                return true;
            }
        }

        // entry not found
        return false;
    }


    //! Return true if name is found in hash table, insert if not found
    bool foundOrInsert(const std::string& name)
    {
        const int hashIdx = hashKeyIndex(name);

        for (hashedEntry* ep = table_[hashIdx]; ep; ep = ep->next_)
        {
            if (name == ep->key_)
            {
                // found - return true
                return true;
            }
        }

        // not found - insert it
        table_[hashIdx] = new hashedEntry(name, table_[hashIdx]);

        // entry not found (but was added) - return false
        return false;
    }

};


/*---------------------------------------------------------------------------*/



#include "wmkdependScanner.h"

namespace wmake {


/*---------------------------------------------------------------------------*\
                           Class Errors Declaration
\*---------------------------------------------------------------------------*/
//! Parser error handing
class Errors
{
public:
	int count;      //!< The number of errors detected

	//! Return a string describing the given error code.
	static std::wstring strerror(int n);

	Errors();               //!< Construct null - start with no errors
	virtual ~Errors();      //!< Destructor
	virtual void clear();   //!< Clear the error count

	//! Handle a general warning 'msg'
	virtual void Warning(const std::wstring& msg);
	//! Handle a general warning 'msg'
	virtual void Warning(int line, int col, const std::wstring& msg);
	//! Handle general error 'msg' (eg, a semantic error)
	virtual void Error(int line, int col, const std::wstring& msg);
	//! Handle syntax error 'n', uses strerror for the message, calls Error()
	virtual void SynErr(int line, int col, int n);
	//! Handle a general exception 'msg'
	virtual void Exception(const std::wstring& msg);

}; // Errors



/*---------------------------------------------------------------------------*\
                           Class Parser Declaration
\*---------------------------------------------------------------------------*/
//! A Coco/R Parser
class Parser
{
	enum {
		_EOF=0,
		_string=1,
		_sqstring=2,
		_package_name=3,
		_package_dir=4,
		maxT = 10    //<! max term (w/o pragmas)
	};
	static const int minErrDist = 2; //!< min. distance before reporting errors

	Token *dummyToken;
	bool deleteErrorsDestruct_; //!< delete the 'errors' member in destructor
	int  errDist;

	void SynErr(int n);         //!< Handle syntax error 'n'
	void Get();
	void Expect(int n);
	bool StartOf(int s);
	void ExpectWeak(int n, int follow);
	bool WeakSeparator(int n, int syFol, int repFol);

public:
	Scanner *scanner;
	Errors  *errors;

	Token *t;                   //!< last recognized token
	Token *la;                  //!< lookahead token

private:

    //! Hash of files already visited
    static StringHashSet visitedFiles_;

    //! Hash of (java) directories already visited
    static StringHashSet visitedDirs_;

    //! Replace all '.' with '/'
    static void dotToSlash(std::string& name);

    //! Import (java) directories
    static void importDir(const std::string& dirName);

    //! Import (java) file
    static void importFile(const std::string& name);

public:
    //! Include directories to search
    static std::list<std::string> includeDirs;

    //! The name of the top-level source file
    static std::string sourceFile;

    //! The name of the top-level dep file
    static std::string depFile;

    //! Add directory to list of visited dirs, thus effectively ignoring it
    static void ignoreDir(const std::string& name);

    //! Include file
    static void includeFile(const std::string& name);

/*---------------------------------------------------------------------------*/

	//! Construct for the specified scanner
	/*!
	 * Use the default error handling, or optionally provide an error
	 * handler, which will not be deleted upon destruction.
	 */
	Parser(Scanner* scan, Errors* err = 0);
	~Parser();
	void Parse();                          //!< Execute the parse operation
	void SemErr(const std::wstring& msg);  //!< Handle semantic error
	bool isUTF8() const;   //!< Return true if scanner buffer is UTF8

	void wmkdepend();

}; // end Parser

} // End namespace

#endif // COCO_wmkdependPARSER_H__

// ************************************************************************* //
