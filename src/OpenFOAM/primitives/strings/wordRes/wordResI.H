/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2021 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

inline const Foam::wordRes& Foam::wordRes::null()
{
    return NullObjectRef<wordRes>();
}


inline Foam::label Foam::wordRes::first_match
(
    const UList<wordRe>& selectors,
    const std::string& text,
    const bool literal
)
{
    label index = 0;
    for (const wordRe& select : selectors)
    {
        if (select.match(text, literal))
        {
            return index;
        }
        ++index;
    }

    return -1;
}


inline Foam::wordRe::compOption Foam::wordRes::found_matched
(
    const UList<wordRe>& selectors,
    const std::string& text
)
{
    auto retval(wordRe::compOption::UNKNOWN);

    for (const wordRe& select : selectors)
    {
        if (select.isLiteral())
        {
            if (select.match(text, true))
            {
                return wordRe::compOption::LITERAL;
            }
        }
        else if
        (
            // Only match regex once
            retval == wordRe::compOption::UNKNOWN
         && select.match(text, false)
        )
        {
            retval = wordRe::compOption::REGEX;
        }
    }

    return retval;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::wordRes::match(const std::string& text, bool literal) const
{
    return (first_match(*this, text, literal) >= 0);
}


inline Foam::wordRe::compOption
Foam::wordRes::matched(const std::string& text) const
{
    return found_matched(*this, text);
}


template<class StringType>
inline Foam::labelList Foam::wordRes::matching
(
    const UList<StringType>& input,
    const bool invert
) const
{
    const label len = input.size();

    labelList indices(len);

    label count = 0;
    for (label i=0; i < len; ++i)
    {
        if (match(input[i]) ? !invert : invert)
        {
            indices[count] = i;
            ++count;
        }
    }
    indices.resize(count);

    return indices;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::wordRes::operator()(const std::string& text) const
{
    return (wordRes::first_match(*this, text) >= 0);
}


inline bool Foam::wordRes::matcher::operator()(const std::string& text) const
{
    return (wordRes::first_match(values, text) >= 0);
}


// ************************************************************************* //
