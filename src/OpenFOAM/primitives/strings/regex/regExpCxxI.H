/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017-2019 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

inline std::regex::flag_type Foam::regExpCxx::syntax()
{
    // 0 = extended, 1 = ECMAScript
    return
    (
        regExpCxx::grammar
      ? std::regex::ECMAScript
      : std::regex::extended
    );
}


inline bool Foam::regExpCxx::meta(const char c)
{
    return
    (
        (c == '.')                           // any character
     || (c == '*' || c == '+' || c == '?')   // quantifiers
     || (c == '(' || c == ')' || c == '|')   // grouping/branching
     || (c == '[' || c == ']')               // range
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::regExpCxx::regExpCxx()
:
    re_(),
    ok_(false)
{}


inline Foam::regExpCxx::regExpCxx(const regExpCxx& rgx)
:
    re_(rgx.re_),
    ok_(rgx.ok_)
{}


inline Foam::regExpCxx::regExpCxx(regExpCxx&& rgx)
:
    re_(std::move(rgx.re_)),
    ok_(rgx.ok_)
{
    rgx.ok_ = false;
}


inline Foam::regExpCxx::regExpCxx(const char* pattern)
:
    re_(),
    ok_(false)
{
    set(pattern, false);
}


inline Foam::regExpCxx::regExpCxx(const std::string& pattern)
:
    re_(),
    ok_(false)
{
    set(pattern, false);
}


inline Foam::regExpCxx::regExpCxx(const char* pattern, bool ignoreCase)
:
    re_(),
    ok_(false)
{
    set(pattern, ignoreCase);
}


inline Foam::regExpCxx::regExpCxx(const std::string& pattern, bool ignoreCase)
:
    re_(),
    ok_(false)
{
    set(pattern, ignoreCase);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline bool Foam::regExpCxx::empty() const
{
    return !ok_;
}


inline bool Foam::regExpCxx::exists() const
{
    return ok_;
}


inline unsigned Foam::regExpCxx::ngroups() const
{
    return ok_ ? re_.mark_count() : 0;
}


inline bool Foam::regExpCxx::nocase() const
{
    return ok_ && ((re_.flags() & std::regex::icase) == std::regex::icase);
}


inline bool Foam::regExpCxx::clear()
{
    if (ok_)
    {
        re_.assign("");
        ok_ = false;

        return true;
    }

    return false;
}


inline void Foam::regExpCxx::swap(regExpCxx& rgx)
{
    re_.swap(rgx.re_);
    std::swap(ok_, rgx.ok_);
}


inline std::string::size_type
Foam::regExpCxx::find(const std::string& text) const
{
    std::smatch mat;
    if (!text.empty() && std::regex_search(text, mat, re_))
    {
        return mat.position(0);
    }

    return std::string::npos;
}


inline bool Foam::regExpCxx::search(const std::string& text) const
{
    return (ok_ && !text.empty() && std::regex_search(text, re_));
}


inline bool Foam::regExpCxx::match(const std::string& text) const
{
    return (ok_ && !text.empty() && std::regex_match(text, re_));
}


inline bool Foam::regExpCxx::match
(
    const std::string& text,
    std::smatch& matches
) const
{
    return std::regex_match(text, matches, re_);
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

inline bool Foam::regExpCxx::operator()(const std::string& text) const
{
    return match(text);
}


inline void Foam::regExpCxx::operator=(const regExpCxx& rgx)
{
    re_ = rgx.re_;
    ok_ = rgx.ok_;
}


inline void Foam::regExpCxx::operator=(regExpCxx&& rgx)
{
    clear();
    swap(rgx);
}


inline void Foam::regExpCxx::operator=(const char* pattern)
{
    set(pattern);
}


inline void Foam::regExpCxx::operator=(const std::string& pattern)
{
    set(pattern);
}


// ************************************************************************* //
