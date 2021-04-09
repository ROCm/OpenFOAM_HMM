/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2021 OpenCFD Ltd.
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
    Foam::CStringList

Description
    An adapter for copying a list of C++ strings into a list of C-style
    strings for passing to C code that expects argc/argv parameters.

    In addition to providing a C-compatible list of C-strings,
    the string lists are flattened into a single string of data that can be
    also be passed en mass.

    Example use:
    \code
        wordList myStrings;  ...
        CStringList cstr(myStrings);

        // pass as argc, argv:
        someMain(cstr.size(), cstr.strings());

        // access the raw characters:
        os.write(cstr.data(), cstr.length());
    \endcode

\*---------------------------------------------------------------------------*/

#ifndef CStringList_H
#define CStringList_H

#include "fileNameList.H"
#include "stringList.H"
#include "wordList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class StringType> class SubStrings;

/*---------------------------------------------------------------------------*\
                         Class CStringList Declaration
\*---------------------------------------------------------------------------*/

class CStringList
{
    // Private Data

        //- Number of strings
        int argc_;

        //- Overall length of the raw content
        //  Does not include the final nul-character
        size_t len_;

        //- List of strings, including trailing nullptr
        char** argv_;

        //- Flattened content with interspersed nul-characters
        char* data_;


    // Private Member Functions

        //- Copy string characters into dest as NUL-terminated string.
        //  Forces conversion of std::sub_match to string form
        //
        //  \return the location one-past the end of dest, which is used
        //      for the next destination
        static inline char* stringCopy(char *dest, const std::string& str);

        //- Copy the input list of strings.
        //  \return number of arguments (argc)
        template<class ListType>
        int resetContent(const ListType& input);


        //- No copy construct
        CStringList(const CStringList&) = delete;

        //- No copy assignment
        void operator=(const CStringList&) = delete;


public:

    // Constructors

        //- Default construct, adding content later (via reset).
        inline CStringList();

        //- Copy construct from a list of strings
        //  Copies the input characters.
        template<class StringType>
        inline explicit CStringList(const UList<StringType>& input);

        //- Copy construct from a list of sub-string references
        //  Copies the input characters.
        template<class StringType>
        inline explicit CStringList(const SubStrings<StringType>& input);


    //- Destructor. Invokes clear() to free memory.
    inline ~CStringList();


    // Public Members

        //- Count the number of parameters until the first nullptr
        //  Return 0 if argv is nullptr.
        static inline int count(const char * const argv[]);


    // Access

        //- True if the size is zero.
        inline bool empty() const noexcept;

        //- Return the number of C-strings (ie, argc)
        inline int size() const noexcept;

        //- Return the list of C-strings (ie, argv)
        //  The position at argc is a nullptr
        inline char** strings() const;

        //- Return the sublist of C-strings (ie, argv) starting at the
        //- specified offset.
        //  \param start the offset, must be less than argc
        inline char** strings(int start) const;


        //- Overall length of the flattened character (data) content
        inline size_t length() const;

        //- The flattened character content, with interspersed nul-chars
        inline const char* data() const;


    // Edit

        //- Clear contents and free memory
        inline void clear();

        //- Copy the input list of strings.
        //  \return number of arguments (argc)
        template<class StringType>
        inline int reset(const UList<StringType>& input);

        //- Copy the input list of strings.
        //  \return number of arguments (argc)
        template<class StringType>
        inline int reset(const SubStrings<StringType>& input);


    // Other

        //- Create a list from argc/argv parameters.
        //  A null pointer for argv is permissible when argc is zero.
        template<class StringType>
        static List<StringType> asList(int argc, const char * const argv[]);

        //- Create a list from a nullptr-terminated list of argv parameters.
        //  Using a nullptr for argv is permissible.
        template<class StringType>
        static inline List<StringType> asList(const char * const argv[]);


    // Member Operators

        //- Return element at the given index. No bounds checking.
        inline const char* operator[](int i) const;
};


// IOstream Operators

//- Output space-separated list
Ostream& operator<<(Ostream& os, const CStringList& list);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CStringListI.H"

#ifdef NoRepository
#   include "CStringListTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
