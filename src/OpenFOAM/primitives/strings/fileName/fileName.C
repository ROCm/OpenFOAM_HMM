/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "fileName.H"
#include "wordList.H"
#include "DynamicList.H"
#include "debug.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::fileName::typeName = "fileName";
int Foam::fileName::debug(debug::debugSwitch(fileName::typeName, 0));
const Foam::fileName Foam::fileName::null;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileName::fileName(const wordList& lst)
{
    forAll(lst, elemI)
    {
        operator=((*this)/lst[elemI]);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName::Type Foam::fileName::type() const
{
    return ::Foam::type(*this);
}


bool Foam::fileName::exists() const
{
    return ::Foam::exists(*this);
}


bool Foam::fileName::isDir() const
{
    return ::Foam::dir(*this);
}


bool Foam::fileName::isFile() const
{
    return ::Foam::file(*this);
}


//  Return file name (part beyond last /)
//
//  behaviour compared to /usr/bin/basename:
//    input           name()          basename
//    -----           ------          --------
//    "foo"           "foo"           "foo"
//    "/foo"          "foo"           "foo"
//    "foo/bar"       "bar"           "bar"
//    "/foo/bar"      "bar"           "bar"
//    "/foo/bar/"     ""              "bar"
//
Foam::word Foam::fileName::name() const
{
    size_type i = rfind('/');

    if (i == npos)
    {
        return *this;
    }
    else
    {
        return substr(i+1, npos);
    }
}


//  Return directory path name (part before last /)
//
//  behaviour compared to /usr/bin/dirname:
//    input           path()          dirname
//    -----           ------          -------
//    "foo"           "."             "."
//    "/foo"          "/"             "foo"
//    "foo/bar"       "foo"           "foo"
//    "/foo/bar"      "/foo"          "/foo"
//    "/foo/bar/"     "/foo/bar/"     "/foo"
//
Foam::fileName Foam::fileName::path() const
{
    size_type i = rfind('/');

    if (i == npos)
    {
        return ".";
    }
    else if (i)
    {
        return substr(0, i);
    }
    else
    {
        return "/";
    }
}


//  Return file name without extension (part before last .)
Foam::fileName Foam::fileName::lessExt() const
{
    size_type i = find_last_of("./");

    if (i == npos || i == 0 || operator[](i) == '/')
    {
        return *this;
    }
    else
    {
        return substr(0, i);
    }
}


//  Return file name extension (part after last .)
Foam::word Foam::fileName::ext() const
{
    size_type i = find_last_of("./");

    if (i == npos || i == 0 || operator[](i) == '/')
    {
        return word::null;
    }
    else
    {
        return substr(i+1, npos);
    }
}


// Return the components of the file name as a wordList
// note that concatenating the components will not necessarily retrieve
// the original input fileName
//
//  behaviour
//    input           components()
//    -----           ------
//    "foo"           1("foo")
//    "/foo"          1("foo")
//    "foo/bar"       2("foo", "bar")
//    "/foo/bar"      2("foo", "bar")
//    "/foo/bar/"     2("foo", "bar")
//
Foam::wordList Foam::fileName::components(const char delimiter) const
{
    DynamicList<word> wrdList(20);

    size_type start=0, end=0;

    while ((end = find(delimiter, start)) != npos)
    {
        // avoid empty element (caused by doubled slashes)
        if (start < end)
        {
            wrdList.append(substr(start, end-start));
        }
        start = end + 1;
    }

    // avoid empty trailing element
    if (start < size())
    {
        wrdList.append(substr(start, npos));
    }

    // transfer to wordList
    return wordList(wrdList.xfer());
}


// Return a component of the file name
Foam::word Foam::fileName::component
(
    const size_type cmpt,
    const char delimiter
) const
{
    return components(delimiter)[cmpt];
}



// Return components following the IOobject requirements
//
//  behaviour
//    input               IOobject(instance, local, name)
//    -----               ------
//    "foo"               ("", "", "foo")
//    "foo/bar"           ("foo", "", "bar")
//    "/XXX"              ERROR - no absolute path
//    "foo/bar/"          ERROR - no name
//    "foo/xxx/bar"       ("foo", "xxx", "bar")
//    "foo/xxx/yyy/bar"   ("foo", "xxx/yyy", "bar")
bool Foam::fileName::IOobjectComponents
(
    fileName& instance,
    fileName& local,
    word& name
)
const
{
    instance.clear();
    local.clear();
    name.clear();

    // called with directory
    if (::Foam::dir(*this))
    {
        std::cerr
            << "fileName::IOobjectComponents() called with directory: "
            << this->c_str() << std::endl;
        std::abort();

        return false;
    }

    size_type first = find('/');

    if (first == 0)
    {
        // called with absolute path
        std::cerr
            << "fileName::IOobjectComponents() called with absolute path: "
            << this->c_str() << std::endl;
        std::abort();

        return false;
    }

    if (first == npos)
    {
        // no '/' found - no instance or local

        // check afterwards
        name.string::operator=(*this);
    }
    else
    {
        instance = substr(0, first);

        size_type last = rfind('/');
        if (last > first)
        {
            // with local
            local = substr(first+1, last-first-1);
        }

        // check afterwards
        name.string::operator=(substr(last+1));
    }


    // check for valid (and stripped) name, regardless of the debug level
    if (name.empty() || string::stripInvalid<word>(name))
    {
        std::cerr
            << "fileName::IOobjectComponents() has invalid word for name: "
            << name.c_str() << "\nwhile processing  "
            << this->c_str() << std::endl;
        std::abort();

        return false;
    }

    return true;
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


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::fileName Foam::operator/(const string& a, const string& b)
{
    if (a.size())           // First string non-null
    {
        if (b.size())       // Second string non-null
        {
            return fileName(a + '/' + b);
        }
        else                // Second string null
        {
            return a;
        }
    }
    else                    // First string null
    {
        if (b.size())       // Second string non-null
        {
            return b;
        }
        else                // Second string null
        {
            return fileName();
        }
    }
}


// ************************************************************************* //
