/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2019 OpenCFD Ltd.
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

#include <algorithm>

// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

inline bool Foam::keyType::valid(char c)
{
    return
    (
        !isspace(c)
     && c != '"'   // string quote
     && c != '\''  // string quote
     && c != '/'   // path separator
     && c != ';'   // end statement
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::keyType::keyType()
:
    word(),
    type_(option::LITERAL)
{}


inline Foam::keyType::keyType(const keyType& s)
:
    word(s, false),
    type_(s.type_)
{}


inline Foam::keyType::keyType(const word& s)
:
    word(s, false),
    type_(option::LITERAL)
{}


inline Foam::keyType::keyType(const string& s)
:
    word(s, false),
    type_(option::REGEX)
{}


inline Foam::keyType::keyType(const char* s)
:
    word(s, false),
    type_(option::LITERAL)
{}


inline Foam::keyType::keyType(const std::string& s, option opt)
:
    word(s, false),
    type_(option(opt & 0x0F))
{}


inline Foam::keyType::keyType(keyType&& s)
:
    word(std::move(static_cast<word&>(s)), false),
    type_(s.type_)
{
    s.type_ = option::LITERAL;
}


inline Foam::keyType::keyType(word&& s)
:
    word(std::move(s), false),
    type_(option::LITERAL)
{}


inline Foam::keyType::keyType(string&& s)
:
    word(std::move(s), false),
    type_(option::REGEX)
{}


inline Foam::keyType::keyType(std::string&& s, option opt)
:
    word(std::move(s), false),
    type_(option(opt & 0x0F))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::keyType::isLiteral() const
{
    return (type_ != option::REGEX);
}


inline bool Foam::keyType::isPattern() const
{
    return (type_ & option::REGEX);
}


inline void Foam::keyType::setType(option opt, bool adjust)
{
    opt = option(opt & 0x0F);

    if (type_ != opt)
    {
        // Only strip when debug is active (potentially costly operation)
        if (isPattern() && adjust && word::debug)
        {
            string::stripInvalid<word>(*this);
        }

        type_ = opt;
    }
}


inline bool Foam::keyType::compile()
{
    type_ = option::REGEX;
    return true;
}


inline void Foam::keyType::uncompile(bool adjust)
{
    setType(option::LITERAL, adjust);
}


inline void Foam::keyType::clear()
{
    word::clear();
    type_ = option::LITERAL;
}


inline void Foam::keyType::swap(keyType& s)
{
    if (this == &s)
    {
        return;  // Self-swap is a no-op
    }

    word::swap(static_cast<word&>(s));
    std::swap(type_, s.type_);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline bool Foam::keyType::operator()(const std::string& text) const
{
    return match(text);  // Use smart match
}


inline void Foam::keyType::operator=(const keyType& s)
{
    if (this == &s)
    {
        return;  // Self-assignment is a no-op
    }

    assign(s); // Bypasses char checking
    type_ = s.type_;
}


inline void Foam::keyType::operator=(keyType&& s)
{
    if (this == &s)
    {
        return;  // Self-assignment is a no-op
    }

    clear();
    swap(s);
}


inline void Foam::keyType::operator=(const word& s)
{
    assign(s); // Bypasses char checking
    type_ = option::LITERAL;
}


inline void Foam::keyType::operator=(const string& s)
{
    assign(s); // Bypasses char checking
    type_ = option::REGEX;
}


inline void Foam::keyType::operator=(const char* s)
{
    assign(s); // Bypasses char checking
    type_ = option::LITERAL;
}


// ************************************************************************* //
