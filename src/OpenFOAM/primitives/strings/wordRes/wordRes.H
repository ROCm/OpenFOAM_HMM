/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2018 OpenCFD Ltd.
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
    Foam::wordRes

Description
    A List of wordRe with additional matching capabilities.

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
:
    public List<wordRe>
{
    // Private Methods

        //- Check for any match of text in list of matchers
        inline static bool found_match
        (
            const UList<wordRe>& patterns,
            const std::string& text,
            bool literal=false
        );


public:

    // Public Classes

        //- Functor wrapper for matching against a List of wordRe
        struct matcher
        {
            const UList<wordRe>& values;

            matcher(const UList<wordRe>& list)
            :
                values(list)
            {}

            //- Return true if string matches ANY of the regular expressions
            //  Allows use as a predicate.
            bool operator()(const std::string& text) const
            {
                return found_match(values, text);
            }
        };


    // Factory Methods

        //- Return a null wordRes - a reference to the NullObject
        inline static const wordRes& null();

        //- Return a wordRes with duplicate words filtered out.
        //  No filtering attempted on regular expressions.
        static wordRes uniq(const UList<wordRe>& input);


    // Constructors

        //- Inherit constructors from List of wordRe
        using List<wordRe>::List;


    //- Destructor
    ~wordRes() = default;


    // Member Functions

        //- Filter out duplicate words (inplace).
        //  No filtering attempted on regular expressions.
        void uniq();

        //- Return true if string matches ANY of the regular expressions
        //  Smart match as regular expression or as a string.
        //  Optionally force a literal match only
        inline bool match(const std::string& text, bool literal=false) const;


    // Member operators

        //- Perform smart match on text, as per match()
        //  Allows use as a predicate.
        inline bool operator()(const std::string& text) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "wordResI.H"

#endif

// ************************************************************************* //
