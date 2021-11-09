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

template<class StringType>
inline Foam::word Foam::IOobject::groupName
(
    StringType base,
    const word& group
)
{
    if (group.empty())
    {
        return base;
    }

    return base + ('.' + group);
}


inline Foam::word Foam::IOobject::scopedName
(
    const std::string& scope,
    const word& name
)
{
    if (scope.empty())
    {
        return name;
    }

    return scope + (IOobject::scopeSeparator + name);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// General access

inline const Foam::word& Foam::IOobject::name() const noexcept
{
    return name_;
}


inline Foam::word Foam::IOobject::group() const
{
    return name_.ext();
}


inline Foam::word Foam::IOobject::member() const
{
    return name_.lessExt();
}


inline const Foam::word& Foam::IOobject::headerClassName() const noexcept
{
    return headerClassName_;
}


inline Foam::word& Foam::IOobject::headerClassName() noexcept
{
    return headerClassName_;
}


inline const Foam::string& Foam::IOobject::note() const noexcept
{
    return note_;
}


inline Foam::string& Foam::IOobject::note() noexcept
{
    return note_;
}


inline bool Foam::IOobject::registerObject() const noexcept
{
    return registerObject_;
}


inline bool Foam::IOobject::registerObject(bool on) noexcept
{
    bool old(registerObject_);
    registerObject_ = on;
    return old;
}


inline bool Foam::IOobject::globalObject() const noexcept
{
    return globalObject_;
}


inline bool Foam::IOobject::globalObject(bool on) noexcept
{
    bool old(globalObject_);
    globalObject_ = on;
    return old;
}


inline unsigned Foam::IOobject::labelByteSize() const noexcept
{
    return static_cast<unsigned>(sizeofLabel_);
}


inline unsigned Foam::IOobject::scalarByteSize() const noexcept
{
    return static_cast<unsigned>(sizeofScalar_);
}


// Checks

inline bool Foam::IOobject::isHeaderClassName(const word& clsName) const
{
    return (clsName == headerClassName_);
}


template<class Type>
inline bool Foam::IOobject::isHeaderClassName() const
{
    return (Type::typeName == headerClassName_);
}


// Read/write options

inline Foam::IOobject::readOption Foam::IOobject::readOpt() const noexcept
{
    return rOpt_;
}


inline Foam::IOobject::readOption
Foam::IOobject::readOpt(readOption opt) noexcept
{
    readOption old(rOpt_);
    rOpt_ = opt;
    return old;
}


inline Foam::IOobject::writeOption Foam::IOobject::writeOpt() const noexcept
{
    return wOpt_;
}


inline Foam::IOobject::writeOption
Foam::IOobject::writeOpt(writeOption opt) noexcept
{
    writeOption old(wOpt_);
    wOpt_ = opt;
    return old;
}


// Path components

inline const Foam::fileName& Foam::IOobject::instance() const noexcept
{
    return instance_;
}


inline Foam::fileName& Foam::IOobject::instance() noexcept
{
    return instance_;
}


inline const Foam::fileName& Foam::IOobject::local() const noexcept
{
    return local_;
}


inline Foam::fileName Foam::IOobject::objectPath() const
{
    return path()/name();
}


// Error Handling

inline bool Foam::IOobject::good() const noexcept
{
    return objState_ == GOOD;
}


inline bool Foam::IOobject::bad() const noexcept
{
    return objState_ == BAD;
}


// ************************************************************************* //
