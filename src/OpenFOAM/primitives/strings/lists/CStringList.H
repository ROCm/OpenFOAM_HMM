/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
     \\/     M anipulation  |
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

/*---------------------------------------------------------------------------*\
                         Class CStringList Declaration
\*---------------------------------------------------------------------------*/

class CStringList
{
    // Private data

        //- Number of strings
        int argc_;

        //- Overall length of the raw content
        //  Does not include the final nul-character
        size_t len_;

        //- List of strings, including trailing NULL pointer
        char** argv_;

        //- Flattened content with interspersed nul-characters
        char* data_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        CStringList(const CStringList&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const CStringList&) = delete;


public:

    // Constructors

        //- Construct empty, adding content later (via reset).
        inline CStringList();


        //- Construct from a list of strings
        //  Copies the input characters.
        template<class StringType>
        CStringList(const UList<StringType>& input);


    //- Destructor
    inline ~CStringList();


    // Public Members

        //- Count the number of parameters until the first NULL pointer.
        //  Return 0 if argv is NULL.
        static inline int count(const char * const argv[]);


    // Access

        //- Return the number of C-strings (ie, argc)
        inline int size() const;

        //- Return the list of C-strings (ie, argv)
        //  The position at argc is a NULL pointer
        inline char** strings() const;


        //- Overall length of the flattened character (data) content
        inline size_t length() const;

        //- The flattened character content, with interspersed nul-chars
        inline char* data() const;


    // Edit

        //- Clear contents and free memory
        inline void clear();

        //- Copy the input list of strings.
        template<class StringType>
        void reset(const UList<StringType>& input);


    // Other

        //- Create a list from argc/argv parameters.
        //  A null pointer for argv is permissible when argc is zero.
        template<class StringType>
        static List<StringType> asList(int argc, const char * const argv[]);

        //- Create a list from a NULL-terminated list of argv parameters.
        //  A null pointer for argv is permissible.
        template<class StringType>
        static inline List<StringType> asList(const char * const argv[]);

};


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
