/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "fileName.H"
#include "wordList.H"
#include "DynamicList.H"
#include "OSspecific.H"
#include "wordRe.H"
#include "fileOperation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::fileName::typeName = "fileName";
int Foam::fileName::debug(debug::debugSwitch(fileName::typeName, 0));
const Foam::fileName Foam::fileName::null;


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::fileName Foam::fileName::validate
(
    const std::string& s,
    const bool doClean
)
{
    fileName out;
    out.resize(s.size());

    char prev = 0;
    std::string::size_type count = 0;

    // Largely as per stripInvalid
    for (auto iter = s.cbegin(); iter != s.cend(); ++iter)
    {
        const char c = *iter;

        if (fileName::valid(c))
        {
            if (doClean && prev == '/' && c == '/')
            {
                // Avoid repeated '/';
                continue;
            }

            // Only track valid chars
            out[count++] = prev = c;
        }
    }

    if (doClean && prev == '/' && count > 1)
    {
        // Avoid trailing '/'
        --count;
    }

    out.resize(count);

    return out;
}


bool Foam::fileName::equals(const std::string& s1, const std::string& s2)
{
    // Do not use (s1 == s2) or s1.compare(s2) first since this would
    // potentially be doing the comparison twice.

    std::string::size_type i1 = 0;
    std::string::size_type i2 = 0;

    const auto n1 = s1.size();
    const auto n2 = s2.size();

    //Info<< "compare " << s1 << " == " << s2 << endl;
    while (i1 < n1 && i2 < n2)
    {
        //Info<< "check '" << s1[i1] << "' vs '" << s2[i2] << "'" << endl;

        if (s1[i1] != s2[i2])
        {
            return false;
        }

        // Increment to next positions and also skip repeated slashes
        do
        {
            ++i1;
        }
        while (s1[i1] == '/');

        do
        {
            ++i2;
        }
        while (s2[i2] == '/');
    }
    //Info<< "return: " << Switch(i1 == n1 && i2 == n2) << endl;

    // Equal if it made it all the way through both strings
    return (i1 == n1 && i2 == n2);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileName::fileName(const UList<word>& lst)
{
    // Estimate overall size
    size_type sz = lst.size();  // Approx number of '/' needed
    for (const word& item : lst)
    {
        sz += item.size();
    }
    reserve(sz);

    sz = 0;
    for (const word& item : lst)
    {
        if (item.size())
        {
            if (sz++) operator+=('/');
            operator+=(item);
        }
    }
}


Foam::fileName::fileName(std::initializer_list<word> lst)
{
    // Estimate overall size
    size_type sz = lst.size();  // Approx number of '/' needed
    for (const word& item : lst)
    {
        sz += item.size();
    }
    reserve(sz);

    sz = 0;
    for (const word& item : lst)
    {
        if (item.size())
        {
            if (sz++) operator+=('/');
            operator+=(item);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName::Type Foam::fileName::type(const bool followLink) const
{
    return ::Foam::type(*this, followLink);
}


Foam::fileName& Foam::fileName::toAbsolute()
{
    fileName& f = *this;

    if (!f.isAbsolute())
    {
        f = cwd()/f;
        f.clean();
    }

    return f;
}


bool Foam::fileName::clean(std::string& str)
{
    // Start with the top slash found - we are never allowed to go above it
    char prev = '/';
    auto top = str.find(prev);

    // No slashes - nothing to do
    if (top == std::string::npos)
    {
        return false;
    }

    // Number of output characters
    std::string::size_type nChar = top+1;

    const string::size_type maxLen = str.size();

    for (string::size_type src = nChar; src < maxLen; /*nil*/)
    {
        const char c = str[src++];

        if (prev == '/')
        {
            // Repeated '/' - skip it
            if (c == '/')
            {
                continue;
            }

            // Could be "/./", "/../" or a trailing "/."
            if (c == '.')
            {
                // Trailing "/." - skip it
                if (src >= maxLen)
                {
                    break;
                }

                // Peek at the next character
                const char c1 = str[src];

                // Found "/./" - skip it
                if (c1 == '/')
                {
                    ++src;
                    continue;
                }

                // Trailing "/.." or intermediate "/../"
                if (c1 == '.' && (src+1 >= maxLen || str[src+1] == '/'))
                {
                    string::size_type parent;

                    // Backtrack to find the parent directory
                    // Minimum of 3 characters:  '/x/../'
                    // Strip it, provided it is above the top point
                    if
                    (
                        nChar > 2
                     && (parent = str.rfind('/', nChar-2)) != string::npos
                     && parent >= top
                    )
                    {
                        nChar = parent + 1;   // Retain '/' from the parent
                        src += 2;
                        continue;
                    }

                    // Bad resolution, eg 'abc/../../'
                    // Retain the sequence, but move the top to avoid it being
                    // considered a valid parent later
                    top = nChar + 2;
                }
            }
        }
        str[nChar++] = prev = c;
    }

    // Remove trailing slash
    if (nChar > 1 && str[nChar-1] == '/')
    {
        nChar--;
    }

    str.resize(nChar);

    return (nChar != maxLen);
}


bool Foam::fileName::clean()
{
    return fileName::clean(*this);
}


Foam::fileName Foam::fileName::clean() const
{
    fileName cleaned(*this);
    fileName::clean(cleaned);
    return cleaned;
}


std::string Foam::fileName::name(const std::string& str)
{
    const auto beg = str.rfind('/');

    if (beg == npos)
    {
        return str;
    }

    return str.substr(beg+1);
}


Foam::word Foam::fileName::name() const
{
    return fileName::name(*this);
}


std::string Foam::fileName::nameLessExt(const std::string& str)
{
    size_type beg = str.rfind('/');
    size_type dot = str.rfind('.');

    if (beg == npos)
    {
        beg = 0;
    }
    else
    {
        ++beg;
    }

    if (dot != npos && dot <= beg)
    {
        dot = npos;
    }

    if (dot == npos)
    {
        return str.substr(beg);
    }

    return str.substr(beg, dot - beg);
}


Foam::word Foam::fileName::nameLessExt() const
{
    return nameLessExt(*this);
}


std::string Foam::fileName::path(const std::string& str)
{
    const auto i = str.rfind('/');

    if (i == npos)
    {
        return ".";
    }
    else if (i)
    {
        return str.substr(0, i);
    }

    return "/";
}


Foam::fileName Foam::fileName::path() const
{
    return path(*this);
}


Foam::fileName Foam::fileName::lessExt() const
{
    const auto i = find_ext();

    if (i == npos)
    {
        return *this;
    }

    return substr(0, i);
}


Foam::word Foam::fileName::ext() const
{
    return string::ext();
}


Foam::fileName& Foam::fileName::ext(const word& ending)
{
    string::ext(ending);
    return *this;
}


bool Foam::fileName::hasExt(const word& ending) const
{
    return string::hasExt(ending);
}


bool Foam::fileName::hasExt(const wordRe& ending) const
{
    return string::hasExt(ending);
}


Foam::wordList Foam::fileName::components(const char delimiter) const
{
    DynamicList<word> wrdList(20);

    size_type beg=0, end=0;

    while ((end = find(delimiter, beg)) != npos)
    {
        // Avoid empty element (caused by doubled slashes)
        if (beg < end)
        {
            wrdList.append(substr(beg, end-beg));
        }
        beg = end + 1;
    }

    // Avoid empty trailing element
    if (beg < size())
    {
        wrdList.append(substr(beg));
    }

    // Transfer to wordList
    return wordList(wrdList.xfer());
}


Foam::word Foam::fileName::component
(
    const size_type cmpt,
    const char delimiter
) const
{
    return components(delimiter)[cmpt];
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::fileName::operator=(const fileName& str)
{
    string::operator=(str);
}


void Foam::fileName::operator=(const word& str)
{
    string::operator=(str);
}


void Foam::fileName::operator=(const string& str)
{
    string::operator=(str);
    stripInvalid();
}


void Foam::fileName::operator=(const std::string& str)
{
    string::operator=(str);
    stripInvalid();
}


void Foam::fileName::operator=(const char* str)
{
    string::operator=(str);
    stripInvalid();
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

Foam::fileName Foam::operator/(const string& a, const string& b)
{
    if (a.size())
    {
        if (b.size())
        {
            // Two non-empty strings: can concatenate
            return fileName(a + '/' + b);
        }

        return a;
    }

    // Or, if the first string is empty

    if (b.size())
    {
        return b;
    }

    // Both strings are empty
    return fileName();
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::search(const word& file, const fileName& directory)
{
    // Search the current directory for the file
    fileNameList files(fileHandler().readDir(directory));
    for (const fileName& item : files)
    {
        if (item == file)
        {
            return directory/item;
        }
    }

    // If not found search each of the sub-directories
    fileNameList dirs(fileHandler().readDir(directory, fileName::DIRECTORY));
    for (const fileName& item : dirs)
    {
        fileName path = search(file, directory/item);
        if (path != fileName::null)
        {
            return path;
        }
    }

    return fileName::null;
}


// ************************************************************************* //
