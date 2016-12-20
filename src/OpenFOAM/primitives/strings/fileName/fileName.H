/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::fileName

Description
    A class for handling file names.

    A fileName is a string of characters without whitespace or quotes.
    A fileName can be
      - constructed from a char*, a string or a word
      - concatenated by adding a '/' separator
      - decomposed into the path, name or component list
      - interrogated for type and access mode

    The string::expand() method expands environment variables, etc,

SourceFiles
    fileName.C
    fileNameIO.C

\*---------------------------------------------------------------------------*/

#ifndef fileName_H
#define fileName_H

#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class T> class List;
template<class T> class UList;
typedef List<word> wordList;

// Forward declaration of friend functions and operators
class wordRe;
class fileName;

Istream& operator>>(Istream&, fileName&);
Ostream& operator<<(Ostream&, const fileName&);


/*---------------------------------------------------------------------------*\
                          Class fileName Declaration
\*---------------------------------------------------------------------------*/

class fileName
:
    public string
{
    // Private Member Functions

        //- Find position of the file extension dot, return npos on failure.
        //  A wrapped version of find_last_of("./") with additional logic.
        inline size_type find_ext() const;

        //- Strip invalid characters
        inline void stripInvalid();


public:

    //- Enumerations to handle file types and modes.
    enum Type
    {
        UNDEFINED,
        FILE,
        DIRECTORY,
        LINK
    };


    // Static data members

        static const char* const typeName;
        static int debug;

        //- An empty fileName
        static const fileName null;


    // Constructors

        //- Construct null
        inline fileName();

        //- Construct as copy
        inline fileName(const fileName& fn);

        //- Construct as copy of word
        inline fileName(const word& s);

        //- Construct as copy of string
        inline fileName(const string& s, const bool doStripInvalid=true);

        //- Construct as copy of std::string
        inline fileName(const std::string& s, const bool doStripInvalid=true);

        //- Construct as copy of character array
        inline fileName(const char* s, const bool doStripInvalid=true);

        //- Construct by concatenating elements of wordList separated by '/'
        explicit fileName(const UList<word>& lst);

        //- Construct by concatenating words separated by '/'
        explicit fileName(std::initializer_list<word> lst);


        //- Construct from Istream
        fileName(Istream& is);


    // Member functions

        //- Is this character valid for a fileName?
        inline static bool valid(char c);

        //- Cleanup file name
        //
        // * Removes repeated slashes
        //       /abc////def        -->   /abc/def
        //
        // * Removes '/./'
        //       /abc/def/./ghi/.   -->   /abc/def/./ghi
        //       abc/def/./         -->   abc/def
        //
        // * Removes '/../'
        //       /abc/def/../ghi/jkl/nmo/..   -->   /abc/ghi/jkl
        //       abc/../def/ghi/../jkl        -->   abc/../def/jkl
        //
        // * Removes trailing '/'
        //
        bool clean();

        //- Cleanup file name
        //  eg, remove repeated slashes, etc.
        fileName clean() const;


        // Interrogation

            //- Return the file type: FILE, DIRECTORY, UNDEFINED or
            //  LINK (only if followLink=false)
            Type type(const bool followLink = true) const;

            //- Return true if file name is absolute
            bool isAbsolute() const;

            //- Convert from relative to absolute
            fileName& toAbsolute();


        // Decomposition

            //- Return file name (part beyond last /)
            //
            //  Behaviour compared to /usr/bin/basename:
            //    Input           name()          basename
            //    -----           ------          --------
            //    "foo"           "foo"           "foo"
            //    "/foo"          "foo"           "foo"
            //    "foo/bar"       "bar"           "bar"
            //    "/foo/bar"      "bar"           "bar"
            //    "/foo/bar/"     ""              "bar"
            //
            word name() const;

            //- Return file name (part beyond last /), subsitute for FOAM_CASE
            string caseName() const;

            //- Return file name, optionally without extension
            word name(const bool noExt) const;

            //- Return directory path name (part before last /)
            //
            //  Behaviour compared to /usr/bin/dirname:
            //    input           path()          dirname
            //    -----           ------          -------
            //    "foo"           "."             "."
            //    "/foo"          "/"             "foo"
            //    "foo/bar"       "foo"           "foo"
            //    "/foo/bar"      "/foo"          "/foo"
            //    "/foo/bar/"     "/foo/bar/"     "/foo"
            //
            fileName path() const;

            //- Return file name without extension (part before last .)
            fileName lessExt() const;

            //- Return file name extension (part after last .)
            word ext() const;

            //- Append a '.' and the ending, and return the object.
            //  The '.' and ending will not be added when the ending is empty,
            //  or when the file name is empty or ended with a '/'.
            fileName& ext(const word& ending);

            //- Return true if it has an extension or simply ends with a '.'
            bool hasExt() const;

            //- Return true if the extension is the same as the given ending.
            bool hasExt(const word& ending) const;

            //- Return true if the extension matches the given ending.
            bool hasExt(const wordRe& ending) const;

            //- Remove extension, returning true if string changed.
            bool removeExt();


            //- Return path components as wordList
            //
            //  Behaviour:
            //    Input           components()
            //    -----           ------
            //    "foo"           1("foo")
            //    "/foo"          1("foo")
            //    "foo/bar"       2("foo", "bar")
            //    "/foo/bar"      2("foo", "bar")
            //    "/foo/bar/"     2("foo", "bar")
            wordList components(const char delimiter = '/') const;

            //- Return a single component of the path
            word component
            (
                const size_type cmpt,
                const char delimiter = '/'
            ) const;


    // Member operators

      // Assignment

        //- Copy, no character validation required
        void operator=(const fileName& str);

        //- Copy, no character validation required
        void operator=(const word& str);

        //- Copy, stripping invalid characters
        void operator=(const string& str);

        //- Copy, stripping invalid characters
        void operator=(const std::string& str);

        //- Copy, stripping invalid characters
        void operator=(const char* str);


    // IOstream operators

        friend Istream& operator>>(Istream& is, fileName& fn);
        friend Ostream& operator<<(Ostream& os, const fileName& fn);
};


//- Assemble words and fileNames as pathnames by adding a '/' separator.
//  No '/' separator is added if either argument is an empty string.
fileName operator/(const string& a, const string& b);


//- Recursively search the given directory for the file
//  returning the path relative to the directory or
//  fileName::null if not found
fileName search(const word& file, const fileName& directory);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fileNameI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
