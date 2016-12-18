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
    Foam::string

Description
    A class for handling character strings derived from std::string.

    Strings may contain any characters and therefore are delimited by quotes
    for IO : "any list of characters".

    Used as a base class for word and fileName.

See also
    Foam::findEtcFile() for information about the site/user OpenFOAM
    configuration directory

SourceFiles
    string.C
    stringIO.C

\*---------------------------------------------------------------------------*/

#ifndef string_H
#define string_H

#include "char.H"
#include "Hasher.H"

#include <string>
#include <cstring>
#include <cstdlib>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class string;
Istream& operator>>(Istream&, string&);
Ostream& operator<<(Ostream&, const string&);
Ostream& operator<<(Ostream&, const std::string&);


/*---------------------------------------------------------------------------*\
                           Class string Declaration
\*---------------------------------------------------------------------------*/

class string
:
    public std::string
{
public:

    // Static data members

        static const char* const typeName;
        static int debug;

        //- An empty string
        static const string null;


    //- Hashing function class, shared by all the derived classes
    class hash
    {
    public:
        hash()
        {}

        inline unsigned operator()(const string&, unsigned seed = 0) const;
    };


    // Constructors

        //- Construct null
        inline string();

        //- Construct from std::string
        inline string(const std::string& str);

        //- Construct as copy of character array
        inline string(const char* str);

        //- Construct as copy of specified number of characters
        inline string(const char* str, const size_type len);

        //- Construct from a single character
        inline string(const char c);

        //- Construct from Istream
        string(Istream& is);


    // Member Functions

        //- Count and return the number of a given character in the string
        size_type count(const char c) const;

        //- Is this string type valid?
        template<class String>
        static inline bool valid(const string& str);

        //- Does this string have particular meta-characters?
        //  The meta characters can be optionally quoted.
        template<class String>
        static inline bool meta(const string& str, const char quote = '\\');

        //- Strip invalid characters from the given string
        template<class String>
        static inline bool stripInvalid(string& str);

        //- Return a valid String from the given string
        template<class String>
        static inline String validate(const string& str);

        //- Return a String with quoted meta-characters from the given string
        template<class String>
        static inline string quotemeta
        (
            const string& str,
            const char quote = '\\'
        );

        //- True when strings match literally
        inline bool match(const std::string& str) const;

        //- Avoid masking the normal std::string replace
        using std::string::replace;

        //- Replace first occurence of sub-string oldStr with newStr
        //  starting at start
        string& replace
        (
            const string& oldStr,
            const string& newStr,
            size_type start = 0
        );

        //- Replace all occurences of sub-string oldStr with newStr
        //  starting at start
        string& replaceAll
        (
            const string& oldStr,
            const string& newStr,
            size_type start = 0
        );

        //- Expand initial tildes and all occurences of environment variables
        //  Expansion includes:
        //  -# environment variables
        //    - "$VAR", "${VAR}"
        //  -# current directory
        //    - leading "./" : the current directory
        //  -# tilde expansion
        //    - leading "~/" : home directory
        //    - leading "~user" : home directory for specified user
        //    - leading "~OpenFOAM" : site/user OpenFOAM configuration directory
        //
        //  Any unknown entries are removed silently if allowEmpty is true
        //  \sa
        //  Foam::findEtcFile
        string& expand(const bool allowEmpty = false);

        //- Remove repeated characters returning true if string changed
        bool removeRepeated(const char character);

        //- Return string with repeated characters removed
        string removeRepeated(const char character) const;

        //- Remove trailing character returning true if string changed
        bool removeTrailing(const char character);

        //- Return string with trailing character removed
        string removeTrailing(const char character) const;

        //- True if the string starts with the given text.
        //  Always true if the given text is empty or if the string
        //  is identical to the given text.
        bool startsWith(const std::string& text) const;

        //- True if the string ends with the given text.
        //  Always true if the given text is empty or if the string
        //  is identical to the given text.
        bool endsWith(const std::string& text) const;


    // Member Operators

        //- Return the sub-string from the i-th character for \a n characters
        inline string operator()
        (
            const size_type i,
            const size_type n
        ) const;

        //- Return the sub-string from the first character for \a n characters
        inline string operator()
        (
            const size_type n
        ) const;


    // IOstream Operators

        friend Istream& operator>>(Istream& is, string& s);
        friend Ostream& operator<<(Ostream& os, const string& s);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "stringI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
