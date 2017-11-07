/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

inline bool Foam::wordRe::meta(const char c)
{
    return regExp::meta(c);
}


inline bool Foam::wordRe::valid(char c)
{
    return keyType::valid(c);
}


inline bool Foam::wordRe::isPattern(const std::string& str)
{
    return string::meta<regExp>(str);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::wordRe::wordRe()
:
    word(),
    re_()
{}


inline Foam::wordRe::wordRe(const wordRe& str)
:
    word(str, false),
    re_()
{
    if (str.isPattern())
    {
        compile();
    }
}


inline Foam::wordRe::wordRe(const keyType& str)
:
    word(str, false),
    re_()
{
    if (str.isPattern())
    {
        compile();
    }
}


inline Foam::wordRe::wordRe(const char* str)
:
    word(str, false),
    re_()
{}


inline Foam::wordRe::wordRe(const std::string& str)
:
    word(str, false),
    re_()
{}


inline Foam::wordRe::wordRe(const string& str)
:
    word(str, false),
    re_()
{}


inline Foam::wordRe::wordRe(const word& str)
:
    word(str, false),
    re_()
{}


inline Foam::wordRe::wordRe(const keyType& str, const compOption opt)
:
    word(str, false),
    re_()
{
    if (str.isPattern())
    {
        compile(opt);
    }
}


inline Foam::wordRe::wordRe(const char* str, const compOption opt)
:
    wordRe(str)
{
    compile(opt);
}


inline Foam::wordRe::wordRe(const std::string& str, const compOption opt)
:
    wordRe(str)
{
    compile(opt);
}


inline Foam::wordRe::wordRe(const string& str, const compOption opt)
:
    wordRe(str)
{
    compile(opt);
}


inline Foam::wordRe::wordRe(const word& str, const compOption opt)
:
    wordRe(str)
{
    compile(opt);
}


inline Foam::wordRe::wordRe(wordRe&& str)
:
    word(std::move(static_cast<word&>(str))),
    re_(std::move(str.re_))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::wordRe::isPattern() const
{
    return re_.exists();
}


inline bool Foam::wordRe::compile(const compOption opt) const
{
    if (opt)
    {
        bool comp = false;

        if (opt & wordRe::REGEX)
        {
            comp = true;
        }
        else if (opt & wordRe::DETECT)
        {
            comp = string::meta<regExp>(*this) || !string::valid<word>(*this);
        }
        else if (opt & wordRe::NOCASE)
        {
            comp = true;
        }

        if (comp)
        {
            return re_.set(*this, (opt & wordRe::NOCASE));
        }
    }

    // Fall-through behaviour - not a regex
    re_.clear();
    return false;
}


inline bool Foam::wordRe::compile() const
{
    return re_.set(*this);
}


inline void Foam::wordRe::uncompile(const bool doStripInvalid) const
{
    if (re_.clear() && doStripInvalid)
    {
        // Skip stripping unless debug is active to avoid costly operations
        if (word::debug)
        {
            string::stripInvalid<word>
            (
                const_cast<word&>(static_cast<const word&>(*this))
            );
        }
    }
}


inline void Foam::wordRe::clear()
{
    word::clear();
    re_.clear();
}


inline bool Foam::wordRe::match(const std::string& text, bool literal) const
{
    if (literal || !re_.exists())
    {
        return !compare(text);  // Compare as literal string
    }
    else
    {
        return re_.match(text); // Match as regex
    }
}


inline Foam::string Foam::wordRe::quotemeta() const
{
    return string::quotemeta<regExp>(*this);
}


inline void Foam::wordRe::set(const std::string& str, const compOption opt)
{
    string::operator=(str);
    compile(opt);
}


inline void Foam::wordRe::set(const char* str, const compOption opt)
{
    string::operator=(str);
    compile(opt);
}


inline void Foam::wordRe::swap(wordRe& str)
{
    word::swap(static_cast<word&>(str));
    re_.swap(str.re_);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::wordRe::operator()(const std::string& text) const
{
    return match(text);
}


inline void Foam::wordRe::operator=(const wordRe& str)
{
    string::operator=(str);
    if (str.isPattern())
    {
        compile();
    }
    else
    {
        re_.clear();
    }
}


inline void Foam::wordRe::operator=(const word& str)
{
    string::operator=(str);
    re_.clear();
}


inline void Foam::wordRe::operator=(const keyType& str)
{
    string::operator=(str);
    if (str.isPattern())
    {
        compile();
    }
    else
    {
        re_.clear();
    }
}


inline void Foam::wordRe::operator=(const string& str)
{
    string::operator=(str);
    compile(wordRe::DETECT);  // Auto-detect regex
}


inline void Foam::wordRe::operator=(const std::string& str)
{
    string::operator=(str);
    compile(wordRe::DETECT);  // Auto-detect regex
}


inline void Foam::wordRe::operator=(const char* str)
{
    string::operator=(str);
    compile(wordRe::DETECT);  // Auto-detect regex
}


inline void Foam::wordRe::operator=(wordRe&& str)
{
    clear();
    swap(str);
}


// ************************************************************************* //
