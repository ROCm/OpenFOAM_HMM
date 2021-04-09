/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

inline bool Foam::wordRe::valid(const char c)
{
    // Also accept '{' and '}' (for regex grouping?)
    return (word::valid(c) || c == '{' || c == '}');
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::wordRe::wordRe()
:
    word(),
    re_()
{}


inline Foam::wordRe::wordRe(const wordRe& str)
:
    word(static_cast<const word&>(str))
{
    if (str.isPattern())
    {
        compile();
    }
}


inline Foam::wordRe::wordRe(wordRe&& str)
:
    word(std::move(static_cast<word&>(str))),
    re_(std::move(str.re_))
{}


inline Foam::wordRe::wordRe(const word& str)
:
    word(str)
{}


inline Foam::wordRe::wordRe(word&& str)
:
    word(std::move(str))
{}


inline Foam::wordRe::wordRe(const std::string& str, const compOption opt)
:
    word(str, false)  // No stripping
{
    if (opt != wordRe::LITERAL)
    {
        compile(opt);
    }
}


inline Foam::wordRe::wordRe(const char* str, const compOption opt)
:
    word(str, false)  // No stripping
{
    if (opt != wordRe::LITERAL)
    {
        compile(opt);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::wordRe::isLiteral() const noexcept
{
    return !re_.exists();
}


inline bool Foam::wordRe::isPattern() const noexcept
{
    return re_.exists();
}


inline bool Foam::wordRe::compile(const compOption opt)
{
    if (opt != wordRe::LITERAL)
    {
        bool comp = false;

        if (opt & wordRe::REGEX)
        {
            comp = true;
        }
        else if (opt & wordRe::DETECT)
        {
            comp = regExp::is_meta(*this) || !string::valid<word>(*this);
        }
        else if (opt & wordRe::ICASE)
        {
            comp = true;
        }

        if (comp)
        {
            return re_.set(*this, (opt & wordRe::ICASE));
        }
    }

    // Fall-through behaviour - not a regex
    re_.clear();
    return false;
}


inline bool Foam::wordRe::compile()
{
    return re_.set(*this);
}


inline void Foam::wordRe::uncompile()
{
    re_.clear();
}


inline void Foam::wordRe::uncompile(bool adjust)
{
    // Only strip when debug is active (potentially costly operation)
    if (adjust && isPattern() && word::debug)
    {
        string::stripInvalid<word>(*this);
    }
    re_.clear();
}


inline void Foam::wordRe::clear()
{
    word::clear();
    re_.clear();
}


inline bool Foam::wordRe::match(const std::string& text, bool literal) const
{
    if (!literal && re_.exists())
    {
        return re_.match(text);  // Match as regex
    }

    return !compare(text);  // Compare as literal
}


inline void Foam::wordRe::set(const std::string& str, const compOption opt)
{
    assign(str);
    compile(opt);
}


inline void Foam::wordRe::set(const char* str, const compOption opt)
{
    assign(str);
    compile(opt);
}


inline void Foam::wordRe::swap(wordRe& str)
{
    if (this == &str)
    {
        return;  // Self-swap is a no-op
    }

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
    if (this == &str)
    {
        return;  // Self-assignment is a no-op
    }

    assign(str);
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
    assign(str);
    re_.clear();
}


inline void Foam::wordRe::operator=(const string& str)
{
    assign(str);
    compile(wordRe::DETECT);  // Auto-detect regex
}


inline void Foam::wordRe::operator=(const std::string& str)
{
    assign(str);
    compile(wordRe::DETECT);  // Auto-detect regex
}


inline void Foam::wordRe::operator=(const char* str)
{
    assign(str);
    compile(wordRe::DETECT);  // Auto-detect regex
}


inline void Foam::wordRe::operator=(wordRe&& str)
{
    if (this == &str)
    {
        return;  // Self-assignment is a no-op
    }

    clear();
    swap(str);
}


// ************************************************************************* //
