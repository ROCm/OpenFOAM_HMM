/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "string.H"
#include "stringOps.H"
#include "word.H"
#include "wordRe.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

const char* const Foam::string::typeName = "string";
int Foam::string::debug(Foam::debug::debugSwitch(string::typeName, 0));
const Foam::string Foam::string::null;


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::word Foam::string::ext() const
{
    const size_type i = find_ext();

    if (i == npos)
    {
        return word::null;
    }
    else
    {
        return substr(i+1, npos);
    }
}


bool Foam::string::ext(const Foam::word& ending)
{
    if (!ending.empty() && !empty() && operator[](size()-1) != '/')
    {
        append(1u, '.');
        append(ending);

        return true;
    }

    return false;
}


bool Foam::string::hasExt(const word& ending) const
{
    size_type i = find_ext();
    if (i == npos)
    {
        return false;
    }

    ++i; // Compare *after* the dot
    return
    (
        // Lengths must match
        ((size() - i) == ending.size())
     && !compare(i, npos, ending)
    );
}


bool Foam::string::hasExt(const wordRe& ending) const
{
    const size_type i = find_ext();
    if (i == npos)
    {
        return false;
    }

    const std::string end = substr(i+1, npos);  // Compare *after* the dot
    return ending.match(end);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::string::size_type Foam::string::count(const char c) const
{
    return stringOps::count(*this, c);
}


Foam::string& Foam::string::replace
(
    const string& oldStr,
    const string& newStr,
    const size_type start
)
{
    size_type pos = start;

    if ((pos = find(oldStr, pos)) != npos)
    {
        std::string::replace(pos, oldStr.size(), newStr);
    }

    return *this;
}


Foam::string& Foam::string::replaceAll
(
    const string& oldStr,
    const string& newStr,
    const size_type start
)
{
    const size_type lenOld = oldStr.size();
    const size_type lenNew = newStr.size();

    if (lenOld)
    {
        for
        (
            size_type pos = start;
            (pos = find(oldStr, pos)) != npos;
            pos += lenNew
        )
        {
            std::string::replace(pos, lenOld, newStr);
        }
    }

    return *this;
}


Foam::string& Foam::string::expand(const bool allowEmpty)
{
    stringOps::inplaceExpand(*this, allowEmpty);
    return *this;
}


bool Foam::string::removeRepeated(const char character)
{
    bool changed = false;

    if (character && find(character) != npos)
    {
        string::size_type nChar = 0;
        iterator outIter = begin();

        char prev = 0;

        for (auto iter = cbegin(); iter != cend(); ++iter)
        {
            const char c = *iter;

            if (prev == c && c == character)
            {
                changed = true;
            }
            else
            {
                *outIter = prev = c;
                ++outIter;
                ++nChar;
            }
        }

        resize(nChar);
    }

    return changed;
}


Foam::string Foam::string::removeRepeated(const char character) const
{
    string str(*this);
    str.removeRepeated(character);
    return str;
}


bool Foam::string::removeTrailing(const char character)
{
    const string::size_type nChar = size();
    if (character && nChar > 1 && operator[](nChar-1) == character)
    {
        resize(nChar-1);
        return true;
    }

    return false;
}


Foam::string Foam::string::removeTrailing(const char character) const
{
    string str(*this);
    str.removeTrailing(character);
    return str;
}


bool Foam::string::removeStart(const std::string& text)
{
    const size_type txtLen = text.size();
    if (!txtLen)
    {
        return true;
    }

    const size_type strLen = this->size();
    if (strLen >= txtLen && !compare(0, txtLen, text))
    {
        this->erase(0, txtLen);
        return true;
    }

    return false;
}


bool Foam::string::removeEnd(const std::string& text)
{
    const size_type txtLen = text.size();
    if (!txtLen)
    {
        return true;
    }

    const size_type strLen = this->size();
    if (strLen >= txtLen && !compare(strLen - txtLen, npos, text))
    {
        this->resize(strLen - txtLen);
        return true;
    }

    return false;
}


bool Foam::string::startsWith(const std::string& text) const
{
    const size_type strLen = this->size();
    const size_type txtLen = text.size();

    return
    (
        !txtLen
     || (strLen >= txtLen && !compare(0, txtLen, text))
    );
}


bool Foam::string::endsWith(const std::string& text) const
{
    const size_type strLen = this->size();
    const size_type txtLen = text.size();

    return
    (
        !txtLen
     || (strLen >= txtLen && !compare(strLen - txtLen, npos, text))
    );
}


// ************************************************************************* //
