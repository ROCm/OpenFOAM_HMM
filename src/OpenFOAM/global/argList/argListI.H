/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "argList.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::word& Foam::argList::executable() const
{
    return executable_;
}


inline const Foam::fileName& Foam::argList::rootPath() const
{
    return rootPath_;
}


inline bool Foam::argList::distributed() const
{
    return distributed_;
}


inline const Foam::fileName& Foam::argList::caseName() const
{
    return case_;
}


inline const Foam::fileName& Foam::argList::globalCaseName() const
{
    return globalCase_;
}


inline const Foam::ParRunControl& Foam::argList::parRunControl() const
{
    return parRunControl_;
}


inline Foam::fileName Foam::argList::path() const
{
    return rootPath()/caseName();
}


inline const Foam::stringList& Foam::argList::args() const
{
    return args_;
}


inline Foam::stringList& Foam::argList::args()
{
    return args_;
}


inline Foam::label Foam::argList::size() const
{
    return args_.size();
}


inline const Foam::HashTable<Foam::string>& Foam::argList::options() const
{
    return options_;
}


inline Foam::HashTable<Foam::string>& Foam::argList::options()
{
    return options_;
}


inline bool Foam::argList::optionFound(const word& optionName) const
{
    return options_.found(optionName);
}


inline Foam::IStringStream Foam::argList::optionLookup
(
    const word& optionName
) const
{
    return IStringStream(options_[optionName]);
}


// * * * * * * * * * * * * Template Specializations  * * * * * * * * * * * * //

namespace Foam
{
    //
    // Specializations for argRead
    //

    template<>
    inline Foam::string
    argList::argRead<Foam::string>(const label argIndex) const
    {
        return args_[argIndex];
    }

    template<>
    inline Foam::word
    argList::argRead<Foam::word>(const label argIndex) const
    {
        return args_[argIndex];
    }

    template<>
    inline Foam::fileName
    argList::argRead<Foam::fileName>(const label argIndex) const
    {
        return args_[argIndex];
    }

    template<>
    inline Foam::label
    argList::argRead<Foam::label>(const label argIndex) const
    {
        return Foam::readLabel(args_[argIndex]);
    }

    template<>
    inline Foam::scalar
    argList::argRead<Foam::scalar>(const label argIndex) const
    {
        return Foam::readScalar(args_[argIndex]);
    }

    //
    // Specializations for optionRead
    //

    template<>
    inline Foam::string
    argList::optionRead<Foam::string>(const word& optionName) const
    {
        return options_[optionName];
    }

    template<>
    inline Foam::word
    argList::optionRead<Foam::word>(const word& optionName) const
    {
        return options_[optionName];
    }

    template<>
    inline Foam::fileName
    argList::optionRead<Foam::fileName>(const word& optionName) const
    {
        return options_[optionName];
    }

    template<>
    inline Foam::label
    argList::optionRead<Foam::label>(const word& optionName) const
    {
        return Foam::readLabel(options_[optionName]);
    }

    template<>
    inline Foam::scalar
    argList::optionRead<Foam::scalar>(const word& optionName) const
    {
        return Foam::readScalar(options_[optionName]);
    }

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline T Foam::argList::argRead(const label argIndex) const
{
    ITstream is(Foam::name(argIndex), args_[argIndex]);

    T val;
    is >> val;

    // Could also check is.nRemainingTokens() to detect trailing rubbish
    return val;
}


template<class T>
inline T Foam::argList::optionRead(const word& optionName) const
{
    ITstream is(optionName, options_[optionName]);

    T val;
    is >> val;

    // Could also check is.nRemainingTokens() to detect trailing rubbish
    return val;
}


template<class T>
inline bool Foam::argList::optionReadIfPresent
(
    const word& optionName,
    T& val
) const
{
    if (optionFound(optionName))
    {
        val = optionRead<T>(optionName);
        return true;
    }

    return false;
}


template<class T>
inline bool Foam::argList::optionReadIfPresent
(
    const word& optionName,
    T& val,
    const T& deflt
) const
{
    if (optionReadIfPresent<T>(optionName, val))
    {
        return true;
    }

    val = deflt;
    return false;
}


template<class T>
inline T Foam::argList::optionLookupOrDefault
(
    const word& optionName,
    const T& deflt
) const
{
    if (optionFound(optionName))
    {
        return optionRead<T>(optionName);
    }

    return deflt;
}


template<class T>
inline Foam::List<T> Foam::argList::optionReadList
(
    const word& optionName
) const
{
    List<T> list;

    ITstream is(optionName, options_[optionName]);

    if (is.size() == 1)
    {
        // Single token - treated like list with one entry
        list.setSize(1);

        is >> list[0];
    }
    else
    {
        is >> list;
    }

    return list;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline const Foam::string& Foam::argList::operator[]
(
    const label argIndex
) const
{
    return args_[argIndex];
}


inline const Foam::string& Foam::argList::operator[]
(
    const word& optionName
) const
{
    return options_[optionName];
}


// ************************************************************************* //
