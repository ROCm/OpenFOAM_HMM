/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2019 OpenCFD Ltd.
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

// Forward declarations
template<class T> class List;
template<class T> class UList;
typedef List<word> wordList;

class wordRe;
class fileName;

/*---------------------------------------------------------------------------*\
                          Class fileName Declaration
\*---------------------------------------------------------------------------*/

class fileName
:
    public string
{
public:

    //- Enumerations to handle directory entry types.
    enum Type
    {
        UNDEFINED = 0,  //!< Undefined type
        FILE = 1,       //!< A file
        DIRECTORY = 2,  //!< A directory
        LINK = 4        //!< A symlink
    };


    // Static Data Members

        //- The typeName
        static const char* const typeName;

        //- Debugging
        static int debug;

        //- Allow space character in fileName. To be used with caution.
        static int allowSpaceInFileName;

        //- An empty fileName
        static const fileName null;


    // Constructors

        //- Construct null
        fileName() = default;

        //- Copy construct
        fileName(const fileName&) = default;

        //- Move construct
        fileName(fileName&&) = default;

        //- Copy construct from word
        inline fileName(const word& s);

        //- Move construct from word
        inline fileName(word&& s);

        //- Copy construct from string
        inline fileName(const string& s, bool doStrip=true);

        //- Move construct from string
        inline fileName(string&& s, bool doStrip=true);

        //- Copy construct from std::string
        inline fileName(const std::string& s, bool doStrip=true);

        //- Move construct from std::string
        inline fileName(std::string&& s, bool doStrip=true);

        //- Copy construct from character array
        inline fileName(const char* s, bool doStrip=true);

        //- Construct by concatenating elements of wordList separated by '/'
        explicit fileName(const UList<word>& list);

        //- Construct by concatenating words separated by '/'
        explicit fileName(std::initializer_list<word> list);

        //- Construct from Istream
        explicit fileName(Istream& is);


    // Member Functions

        //- Is this character valid for a fileName?
        inline static bool valid(char c);

        //- Construct fileName with no invalid characters, possibly applying
        //- other transformations such as changing the path separator,
        //- removing duplicate or trailing slashes, etc.
        static fileName validate(const std::string& s, const bool doClean=true);

        //- Join two strings with a path separator ('/' by default).
        //  No separator is added if either argument is an empty string or
        //  if the arguments already had the path separator at the junction.
        //  Invalid characters are \em not stripped (ie, retained).
        static fileName concat
        (
            const std::string& s1,
            const std::string& s2,
            const char delim = '/'
        );

        //- This is a specialized (possibly slower) version of compare()
        //- that ignores duplicate or trailing slashes.
        static bool equals(const std::string& s1, const std::string& s2);

        //- Strip invalid characters
        inline void stripInvalid();

        //- Cleanup filename
        //
        // Removes trailing \c /
        //   \verbatim
        //       /                 -->   /
        //       /abc/             -->   /abc
        //   \endverbatim
        //
        // Removes repeated slashes
        //   \verbatim
        //       /abc////def        -->   /abc/def
        //   \endverbatim
        //
        // Removes \c /./ (current directory)
        //   \verbatim
        //       /abc/def/./ghi/.   -->   /abc/def/ghi
        //       abc/def/./         -->   abc/def
        //       ./abc/             -->   ./abc
        //   \endverbatim
        //
        // Removes \c /../ (parent directory)
        //   \verbatim
        //       /abc/def/../ghi/jkl/nmo/..   -->   /abc/ghi/jkl
        //       abc/../def/ghi/../jkl        -->   abc/../def/jkl
        //   \endverbatim
        //
        // \return True if the content changed
        static bool clean(std::string& str);


        //- Cleanup filename inplace
        //  \return True if any contents changed
        bool clean();

        //- Cleanup filename
        //  \return cleaned copy of fileName
        fileName clean() const;


    // Interrogation

        //- Return the directory entry type:
        //- UNDEFINED, FILE, DIRECTORY (or LINK).
        //
        //  \param followLink when false it will return LINK for a symlink
        //     rather than following it.
        //  \param checkGzip add an additional test for a gzip FILE
        Type type(bool followLink=true, bool checkGzip=false) const;

        //- Return true if string starts with a '/'
        inline static bool isAbsolute(const std::string& str);

        //- Return true if file name is absolute (starts with a '/')
        inline bool isAbsolute() const;

        //- Convert from relative to absolute
        fileName& toAbsolute();

        //- Return true if string ends with "~", ".bak", ".old", ".save"
        static bool isBackup(const std::string& str);

        //- Return true if file name ends with "~", ".bak", ".old", ".save"
        inline bool isBackup() const;


    // Decomposition

        //- Return basename (part beyond last /), including its extension
        //  The result normally corresponds to a Foam::word
        //
        // Behaviour compared to /usr/bin/basename:
        // \verbatim
        //    input           name()          basename
        //    -----           ------          --------
        //    ""              ""              ""
        //    "abc"           "abc"           "abc"
        //    "/"             ""              "/"
        //    "/abc"          "abc"           "abc"
        //    "abc/def"       "def"           "def"
        //    "/abc/def"      "def"           "def"
        //    "/abc/def/"     ""              "def"
        //    "/abc/../def"   "def"           "def"
        // \endverbatim
        inline static std::string name(const std::string& str);

        //- Return basename (part beyond last /), including its extension
        inline word name() const;

        //- Return basename, without extension
        //  The result normally corresponds to a Foam::word
        static std::string nameLessExt(const std::string& str);

        //- Return basename, without extension
        inline word nameLessExt() const;

        //- Deprecated(2017-03) return basename, optionally without extension
        //  \deprecated(2017-03) - use name() or nameLessExt() methods
        //      which describe their behaviour explicitly
        word name(const bool noExt) const
        {
            return noExt ? this->nameLessExt() : this->name();
        }

        //- Return directory path name (part before last /)
        //  The result normally coresponds to a Foam::fileName
        //
        // Behaviour compared to /usr/bin/dirname:
        // \verbatim
        //    input           path()          dirname
        //    -----           ------          -------
        //    ""              "."             "."
        //    "abc"           "."             "."
        //    "/"             "/"             "/"
        //    "/abc"          "/"             "/"
        //    "abc/def"       "abc"           "abc"
        //    "/abc/def"      "/abc"          "/abc"
        //    "/abc/def/"     "/abc/def"      "/abc"
        //    "/abc/../def"   "/abc/.."       "/abc/.."
        // \endverbatim
        inline static std::string path(const std::string& str);

        //- Return directory path name (part before last /)
        inline fileName path() const;

        //- Return true if it contains a '/' character
        inline bool hasPath() const;

        //- Remove leading path, returning true if string changed.
        inline bool removePath();

        //- Return a relative name by stripping off the parent directory
        //- where possible.
        //
        //  \param parent the parent directory
        //  \param caseTag replace the parent with \<case\> for later
        //      use with expand(), or prefix \<case\> if the file name was
        //      not an absolute location
        fileName relative
        (
            const fileName& parent,
            const bool caseTag = false
        ) const;

        //- Return file name without extension (part before last .)
        inline fileName lessExt() const;

        //- Return file name extension (part after last .)
        inline word ext() const;

        //- Append a '.' and the ending, and return the object.
        //  The '.' and ending will not be added when the ending is empty,
        //  or when the file name is empty or ended with a '/'.
        inline fileName& ext(const word& ending);

        //- Return true if it has an extension or simply ends with a '.'
        inline bool hasExt() const;

        //- Return true if the extension is the same as the given ending.
        inline bool hasExt(const word& ending) const;

        //- Return true if the extension matches the given ending.
        bool hasExt(const wordRe& ending) const;

        //- Remove extension, returning true if string changed.
        inline bool removeExt();


        //- Return path components as wordList
        //
        // Behaviour:
        // \verbatim
        //    input           components()
        //    -----           ------------
        //    ""              ()
        //    "."             (".")
        //    "abc"           ("abc")
        //    "/abc"          ("abc")
        //    "abc/def"       ("abc", "def")
        //    "/abc/def"      ("abc", "def")
        //    "/abc/def/"     ("abc", "def")
        // \endverbatim
        wordList components(const char delim = '/') const;

        //- Return a single component of the path
        word component(const size_type cmpt, const char delim = '/') const;


    // Member Operators

    // Assignment

        //- Copy assignment, no character validation required
        //  Self-assignment is a no-op.
        inline fileName& operator=(const fileName& str);

        //- Move assignment, no character validation required
        //  Self-assignment is a no-op.
        inline fileName& operator=(fileName&& str);

        //- Copy assignment, no character validation required
        inline fileName& operator=(const word& str);

        //- Move assignment, no character validation required
        inline fileName& operator=(word&& str);

        //- Copy assignment, stripping invalid characters
        inline fileName& operator=(const string& str);

        //- Move assignment, stripping invalid characters
        inline fileName& operator=(string&& str);

        //- Copy assignment, stripping invalid characters
        inline fileName& operator=(const std::string& str);

        //- Move assignment, stripping invalid characters
        inline fileName& operator=(std::string&& str);

        //- Copy, stripping invalid characters
        inline fileName& operator=(const char* str);


    // Other operators

        //- Append a path element with '/' separator.
        //  No '/' separator is added if this or the argument are empty.
        fileName& operator/=(const string& other);
};


// IOstream Operators

//- Read operator
Istream& operator>>(Istream& is, fileName& val);

//- Write operator
Ostream& operator<<(Ostream& os, const fileName& val);


// Global Operators

//- Assemble words and fileNames as pathnames by adding a '/' separator.
//  No '/' separator is added if either argument is an empty string.
fileName operator/(const string& s1, const string& s2);


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
