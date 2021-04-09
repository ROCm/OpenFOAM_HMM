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

#include <algorithm>

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::regExpPosix::regExpPosix() noexcept
:
    preg_(nullptr)
{}


inline Foam::regExpPosix::regExpPosix(regExpPosix&& rgx) noexcept
:
    preg_(rgx.preg_)
{
    rgx.preg_ = nullptr;
}


inline Foam::regExpPosix::regExpPosix
(
    const char* pattern,
    const bool ignoreCase
)
:
    preg_(nullptr)
{
    set(pattern, ignoreCase);
}


inline Foam::regExpPosix::regExpPosix
(
    const std::string& pattern,
    const bool ignoreCase
)
:
    preg_(nullptr)
{
    set(pattern, ignoreCase);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

inline Foam::regExpPosix::~regExpPosix()
{
    clear();
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline bool Foam::regExpPosix::empty() const noexcept
{
    return !preg_;
}


inline bool Foam::regExpPosix::exists() const noexcept
{
    return preg_;
}


inline unsigned Foam::regExpPosix::ngroups() const
{
    return preg_ ? preg_->re_nsub : 0;
}


inline bool Foam::regExpPosix::search(const std::string& text) const
{
    return std::string::npos != find(text);
}


inline void Foam::regExpPosix::swap(regExpPosix& rgx)
{
    if (this != &rgx)
    {
        // Self-swap is a no-op
        std::swap(preg_, rgx.preg_);
    }
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

inline bool Foam::regExpPosix::operator()(const std::string& text) const
{
    return match(text);
}


inline void Foam::regExpPosix::operator=(regExpPosix&& rgx)
{
    if (this != &rgx)
    {
        // Self-assignment is a no-op
        clear();
        swap(rgx);
    }
}


inline void Foam::regExpPosix::operator=(const char* pattern)
{
    set(pattern);
}


inline void Foam::regExpPosix::operator=(const std::string& pattern)
{
    set(pattern);
}


// ************************************************************************* //
