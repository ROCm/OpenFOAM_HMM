/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class StringType>
inline Foam::word Foam::IOobject::groupName(StringType name, const word& group)
{
    if (group.empty())
    {
        return name;
    }

    return name + ('.' + group);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// General access

inline const Foam::word& Foam::IOobject::name() const
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


inline const Foam::word& Foam::IOobject::headerClassName() const
{
    return headerClassName_;
}


inline Foam::word& Foam::IOobject::headerClassName()
{
    return headerClassName_;
}


inline const Foam::string& Foam::IOobject::note() const
{
    return note_;
}


inline Foam::string& Foam::IOobject::note()
{
    return note_;
}


inline bool Foam::IOobject::registerObject() const
{
    return registerObject_;
}


inline bool& Foam::IOobject::registerObject()
{
    return registerObject_;
}


inline bool Foam::IOobject::globalObject() const
{
    return globalObject_;
}


inline bool& Foam::IOobject::globalObject()
{
    return globalObject_;
}


inline unsigned Foam::IOobject::labelByteSize() const
{
    return labelByteSize_;
}


inline unsigned Foam::IOobject::scalarByteSize() const
{
    return scalarByteSize_;
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

inline Foam::IOobject::readOption Foam::IOobject::readOpt() const
{
    return rOpt_;
}


inline Foam::IOobject::readOption& Foam::IOobject::readOpt()
{
    return rOpt_;
}


inline Foam::IOobject::writeOption Foam::IOobject::writeOpt() const
{
    return wOpt_;
}


inline Foam::IOobject::writeOption& Foam::IOobject::writeOpt()
{
    return wOpt_;
}


// Path components

inline const Foam::fileName& Foam::IOobject::instance() const
{
    return instance_;
}


inline Foam::fileName& Foam::IOobject::instance()
{
    return instance_;
}


inline const Foam::fileName& Foam::IOobject::local() const
{
    return local_;
}


inline Foam::fileName Foam::IOobject::objectPath() const
{
    return path()/name();
}


// Error Handling

inline bool Foam::IOobject::good() const
{
    return objState_ == GOOD;
}


inline bool Foam::IOobject::bad() const
{
    return objState_ == BAD;
}


// ************************************************************************* //
