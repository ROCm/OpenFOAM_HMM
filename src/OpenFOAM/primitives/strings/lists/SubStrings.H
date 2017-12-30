/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::SubStrings

Description
    Sub-ranges of a string with a structure similar to std::match_results,
    but without the underlying regular expression matching.

\*---------------------------------------------------------------------------*/

#ifndef SubStrings_H
#define SubStrings_H

#include <string>
#include <regex>
#include <vector>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class SubStrings Declaration
\*---------------------------------------------------------------------------*/

template<class String>
class SubStrings
:
    public std::vector<std::sub_match<typename String::const_iterator>>
{
public:

    // Typedefs

        //- The element type
        using value_type =
            typename std::sub_match<typename String::const_iterator>;

        //- The const_iterator for the underlying string type
        using string_iterator = typename String::const_iterator;


    // Constructors

        //- Construct null
        SubStrings()
        {}


    // Member Functions

        //- The total length of all sub-elements.
        //  Use size() for the number elements.
        std::string::size_type length() const
        {
            std::string::size_type len = 0;

            for (const auto& elem : *this)
            {
                len += elem.length();
            }

            return len;
        }


        //- Append sub-string defined by begin/end iterators
        void append(string_iterator b, string_iterator e)
        {
            value_type range;
            range.first = b;
            range.second = e;
            range.matched = true;

            this->push_back(range);
        }


        //- Const reference to the first element,
        //  for consistency with other OpenFOAM containers
        auto first() const -> decltype(this->front())
        {
            return this->front();
        }


        //- Const reference to the last element,
        //  for consistency with other OpenFOAM containers
        auto last() const -> decltype(this->back())
        {
            return this->back();
        }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
