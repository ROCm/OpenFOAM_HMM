/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::wordRes

Description
    A wrapper for matching a std::string against a wordRe list.

Note
    The constructor should remain non-explicit. This allows automatic
    conversion from UList\<wordRe\> to wordRes in search functions.

SourceFiles
    wordResI.H
    wordRes.C

\*---------------------------------------------------------------------------*/

#ifndef wordRes_H
#define wordRes_H

#include "wordReList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class wordRes Declaration
\*---------------------------------------------------------------------------*/

class wordRes
{
    // Private data

        //- Reference to underlying list
        const UList<wordRe>& list_;

public:

    // STL type definitions

        //- Type of values the list contains
        typedef wordRe value_type;


    // Constructors

        //- Construct from a list of wordRe
        inline wordRes(const UList<wordRe>& list);


    // Static Constructors, Helpers

        //- Return a wordReList with duplicate words filtered out.
        //  No filtering is done on regular expressions.
        static wordReList uniq(const UList<wordRe>& input);

        //- Inplace subset of wordReList with duplicate words filtered out.
        //  No filtering is done on regular expressions.
        static void inplaceUniq(wordReList& input);


    // Member Functions

      // Access

        //- The number of elements in the list
        inline label size() const;

        //- True if the list is empty
        inline bool empty() const;


      // Searching

        //- Return true if string matches any of the regular expressions
        //  Smart match as regular expression or as a string.
        //  Optionally specify a literal match only.
        inline bool match
        (
            const std::string& text,
            const bool literal = false
        ) const;


    // Member operators

        //- Return underlying list of wordRe
        inline const UList<wordRe>& operator()() const;

        //- Perform smart match on text, as per match()
        inline bool operator()(const std::string& text) const;

        //- Return element of constant list
        inline const wordRe& operator[](const label i) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "wordResI.H"

#endif

// ************************************************************************* //
