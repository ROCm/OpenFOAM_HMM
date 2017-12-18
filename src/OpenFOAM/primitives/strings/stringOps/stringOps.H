/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

Namespace
    Foam::stringOps

Description
    Collection of static functions to do various simple string-related
    operations

SourceFiles
    stringOps.C
    stringOpsTemplates.C

\*---------------------------------------------------------------------------*/
#ifndef stringOps_H
#define stringOps_H

#include "string.H"
#include "SubStrings.H"
#include "word.H"
#include "dictionary.H"
#include "HashTable.H"
#include "stringOpsSort.H"
#include "wordRes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Namespace stringOps Declaration
\*---------------------------------------------------------------------------*/

namespace stringOps
{
    //- Count the number of occurences of the specified character
    std::string::size_type count(const std::string& str, const char c);

    //- Count the number of occurences of the specified character
    //  Correctly handles nullptr.
    std::string::size_type count(const char* str, const char c);

    //- Return true if text matches one of the regular expressions.
    //  Simply forwards a wordReList to a wordRes for the matching.
    inline bool match(const wordReList& patterns, const std::string& text)
    {
        return wordRes(patterns).match(text);
    }


    //- Expand occurences of variables according to the mapping
    //  Expansion includes:
    //  -# variables
    //    - "$VAR", "${VAR}"
    //
    //  Supports default and alternative values as per the POSIX shell.
    //  \code
    //      a)  "${parameter:-defValue}"
    //      b)  "${parameter:+altValue}"
    //  \endcode
    //  a) If parameter is unset or null, the \c defValue is substituted.
    //  Otherwise, the value of parameter is substituted.
    //
    //  b) If parameter is unset or null, nothing is substituted.
    //  Otherwise the \c altValue is substituted.
    //
    //  - Any unknown entries are removed silently.
    //  - Malformed entries (eg, brace mismatch, sigil followed by bad chars)
    //  are left as is.
    //
    //  \note the leading sigil can be changed to avoid conflicts with other
    //  string expansions
    string expand
    (
        const string& original,
        const HashTable<string, word, string::hash>& mapping,
        const char sigil = '$'
    );


    //- Inplace expand occurences of variables according to the mapping
    //  Expansion includes:
    //  -# variables
    //    - "$VAR", "${VAR}"
    //
    //  Supports default and alternative values as per the POSIX shell.
    //  \code
    //      a)  "${parameter:-defValue}"
    //      b)  "${parameter:+altValue}"
    //  \endcode
    //  a) If parameter is unset or null, the \c defValue is substituted.
    //  Otherwise, the value of parameter is substituted.
    //
    //  b) If parameter is unset or null, nothing is substituted.
    //  Otherwise the \c altValue is substituted.
    //
    //  - Any unknown entries are removed silently.
    //  - Malformed entries (eg, brace mismatch, sigil followed by bad chars)
    //  are left as is.
    //
    //  \note the leading sigil can be changed to avoid conflicts with other
    //  string expansions
    void inplaceExpand
    (
        std::string& s,
        const HashTable<string, word, string::hash>& mapping,
        const char sigil = '$'
    );

    //- Expand occurences of variables according to the dictionary
    //  Expansion includes:
    //  -# variables
    //    - "$VAR", "${VAR}"
    //
    //  Any unknown entries are left as-is
    //
    //  \note the leading sigil can be changed to avoid conflicts with other
    //  string expansions
    string expand
    (
        const string& original,
        const dictionary& dict,
        const char sigil = '$'
    );


    //- Get dictionary or (optionally) environment variable
    //
    //  The environment variable lookup supports default and alternative
    //  values as per the POSIX shell.
    //  \code
    //      ${parameter:-defValue}
    //      ${parameter:+altValue}
    //  \endcode
    string getVariable
    (
        const word& name,
        const dictionary& dict,
        const bool allowEnvVars,
        const bool allowEmpty
    );


    //- Recursively expands (dictionary or environment) variable
    //  starting at index in string. Updates index.
    string expand
    (
        const string& s,
        std::string::size_type& index,
        const dictionary& dict,
        const bool allowEnvVars,
        const bool allowEmpty
    );


    //- Inplace expand occurences of variables according to the dictionary
    //  and optionally environment variables
    //  Expansion includes:
    //  -# variables
    //    - "$VAR", "${VAR}"
    //
    //  with the "${}" syntax doing a recursive substitution.
    //  Any unknown entries are left as-is
    //
    //  \note the leading sigil can be changed to avoid conflicts with other
    //  string expansions
    void inplaceExpand
    (
        std::string& s,
        const dictionary& dict,
        const bool allowEnvVars,
        const bool allowEmpty,
        const char sigil = '$'
    );


    //- Inplace expand occurences of variables according to the dictionary
    //  Expansion includes:
    //  -# variables
    //    - "$VAR", "${VAR}"
    //
    //  Any unknown entries are left as-is
    //
    //  \note the leading sigil can be changed to avoid conflicts with other
    //  string expansions
    void inplaceExpand
    (
        std::string& s,
        const dictionary& dict,
        const char sigil = '$'
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
    //  Supports default and alternative values as per the POSIX shell.
    //  \code
    //      a)  "${parameter:-defValue}"
    //      b)  "${parameter:+altValue}"
    //  \endcode
    //  a) If parameter is unset or null, the \c defValue is substituted.
    //  Otherwise, the value of parameter is substituted.
    //
    //  b) If parameter is unset or null, nothing is substituted.
    //  Otherwise the \c altValue is substituted.
    //
    //  - Any unknown entries are removed silently, if allowEmpty is true.
    //  - Malformed entries (eg, brace mismatch, sigil followed by bad chars)
    //  are left as is.
    //
    //  \sa
    //  Foam::findEtcFile
    string expand
    (
        const string& original,
        const bool allowEmpty = false
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
    //  Supports default and alternative values as per the POSIX shell.
    //  \code
    //      a)  "${parameter:-defValue}"
    //      b)  "${parameter:+altValue}"
    //  \endcode
    //  a) If parameter is unset or null, the \c defValue is substituted.
    //  Otherwise, the value of parameter is substituted.
    //
    //  b) If parameter is unset or null, nothing is substituted.
    //  Otherwise the \c altValue is substituted.
    //
    //  - Any unknown entries are removed silently if allowEmpty is true.
    //  - Malformed entries (eg, brace mismatch, sigil followed by bad chars)
    //  are left as is.
    //
    //  \sa
    //  Foam::findEtcFile
    void inplaceExpand
    (
        std::string& s,
        const bool allowEmpty = false
    );


    //- Replace environment variable contents with its name.
    //  This is essentially the inverse operation for inplaceExpand.
    //  Return true if a replacement was successful.
    bool inplaceReplaceVar(std::string& s, const word& varName);


    //- Return string trimmed of leading whitespace
    string trimLeft(const string& s);

    //- Trim leading whitespace inplace
    void inplaceTrimLeft(std::string& s);

    //- Return string trimmed of trailing whitespace
    string trimRight(const string& s);

    //- Trim trailing whitespace inplace
    void inplaceTrimRight(std::string& s);

    //- Return string trimmed of leading and trailing whitespace
    string trim(const string& original);

    //- Trim leading and trailing whitespace inplace
    void inplaceTrim(std::string& s);


    //- Return string transformed with std::tolower on each character
    string lower(const string& original);

    //- Inplace transform string with std::tolower on each character
    void inplaceLower(std::string& s);

    //- Return string transformed with std::toupper on each character
    string upper(const string& original);

    //- Inplace transform string with std::toupper on each character
    void inplaceUpper(std::string& s);


    //- Using printf-formatter for a word representation of the primitive.
    //  The representation is not checked for valid word characters -
    //  it is assumed that the caller knows what they are doing
    template<class PrimitiveType>
    Foam::word name(const char* fmt, const PrimitiveType& val);

    //- Using printf-formatter for a word representation of the primitive.
    //  The representation is not checked for valid word characters -
    //  it is assumed that the caller knows what they are doing
    template<class PrimitiveType>
    Foam::word name(const std::string& fmt, const PrimitiveType& val);


    //- Split string into sub-strings at the delimiter character.
    //  Empty sub-strings are normally suppressed.
    //  Behaviour is ill-defined if delim is a NUL character.
    template<class StringType>
    Foam::SubStrings<StringType> split
    (
        const StringType& str,
        const char delim,
        const bool keepEmpty = false
    );

    //- Split string into sub-strings using delimiter string.
    //  Empty sub-strings are normally suppressed.
    template<class StringType>
    Foam::SubStrings<StringType> split
    (
        const StringType& str,
        const std::string& delim,
        const bool keepEmpty = false
    );

    //- Split string into sub-strings using any characters in delimiter.
    //  Empty sub-strings are normally suppressed.
    //  Behaviour is ill-defined if delim is an empty string.
    template<class StringType>
    Foam::SubStrings<StringType> splitAny
    (
        const StringType& str,
        const std::string& delim
    );

    //- Split string into sub-strings using a fixed field width.
    //  Behaviour is ill-defined if width is zero.
    //  \param str the string to be split
    //  \param width the fixed field width for each sub-string
    //  \param start the optional offset of where to start the splitting.
    //      Any text prior to start is ignored in the operation.
    template<class StringType>
    Foam::SubStrings<StringType> splitFixed
    (
        const StringType& str,
        const std::string::size_type width,
        const std::string::size_type start = 0
    );

    //- Split string into sub-strings at whitespace (TAB, NL, VT, FF, CR, SPC)
    //  Empty sub-strings are suppressed.
    template<class StringType>
    Foam::SubStrings<StringType> splitSpace
    (
        const StringType& str
    );


} // End namespace stringOps

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "stringOpsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
