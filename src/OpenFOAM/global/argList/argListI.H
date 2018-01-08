/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017-2018 OpenCFD Ltd.
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


inline Foam::label Foam::argList::size() const
{
    return args_.size();
}


inline const Foam::stringList& Foam::argList::args() const
{
    return args_;
}


inline Foam::stringList& Foam::argList::args()
{
    return args_;
}


inline const Foam::HashTable<Foam::string>& Foam::argList::options() const
{
    return options_;
}


inline Foam::HashTable<Foam::string>& Foam::argList::options()
{
    return options_;
}


inline bool Foam::argList::found(const word& optName) const
{
    return options_.found(optName);
}


inline Foam::IStringStream Foam::argList::lookup(const word& optName) const
{
    return IStringStream(options_[optName]);
}


// * * * * * * * * * * * * Template Specializations  * * * * * * * * * * * * //

namespace Foam
{
    template<> inline int32_t argList::read<int32_t>(const label index) const
    {
        return Foam::readInt32(args_[index]);
    }

    template<> inline int64_t argList::read<int64_t>(const label index) const
    {
        return Foam::readInt64(args_[index]);
    }

    template<> inline float argList::read<float>(const label index) const
    {
        return Foam::readFloat(args_[index]);
    }

    template<> inline double argList::read<double>(const label index) const
    {
        return Foam::readDouble(args_[index]);
    }


    template<> inline int32_t argList::opt<int32_t>(const word& optName) const
    {
        return Foam::readInt32(options_[optName]);
    }

    template<> inline int64_t argList::opt<int64_t>(const word& optName) const
    {
        return Foam::readInt64(options_[optName]);
    }

    template<> inline float argList::opt<float>(const word& optName) const
    {
        return Foam::readFloat(options_[optName]);
    }

    template<> inline double argList::opt<double>(const word& optName) const
    {
        return Foam::readDouble(options_[optName]);
    }


    template<>
    inline string argList::read<Foam::string>(const label index) const
    {
        return args_[index];
    }

    template<>
    inline word argList::read<Foam::word>(const label index) const
    {
        return args_[index];
    }

    template<>
    inline fileName argList::read<Foam::fileName>(const label index) const
    {
        return args_[index];
    }


    template<>
    inline string argList::opt<Foam::string>(const word& optName) const
    {
        return options_[optName];
    }

    template<>
    inline word argList::opt<Foam::word>(const word& optName) const
    {
        return options_[optName];
    }

    template<>
    inline fileName argList::opt<Foam::fileName>(const word& optName) const
    {
        return options_[optName];
    }

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline T Foam::argList::read(const label index) const
{
    ITstream is(Foam::name(index), args_[index]);

    T val;
    is >> val;

    // Could also check is.nRemainingTokens() to detect trailing rubbish
    return val;
}


template<class T>
inline T Foam::argList::opt(const word& optName) const
{
    ITstream is(optName, options_[optName]);

    T val;
    is >> val;

    // Could also check is.nRemainingTokens() to detect trailing rubbish
    return val;
}


template<class T>
inline bool Foam::argList::readIfPresent
(
    const word& optName,
    T& val
) const
{
    if (found(optName))
    {
        val = opt<T>(optName);
        return true;
    }

    return false;
}


template<class T>
inline bool Foam::argList::readIfPresent
(
    const word& optName,
    T& val,
    const T& deflt
) const
{
    if (readIfPresent<T>(optName, val))
    {
        return true;
    }

    val = deflt;
    return false;
}


template<class T>
inline T Foam::argList::lookupOrDefault
(
    const word& optName,
    const T& deflt
) const
{
    if (found(optName))
    {
        return opt<T>(optName);
    }

    return deflt;
}


template<class T>
inline Foam::List<T> Foam::argList::readList(const word& optName) const
{
    List<T> list;

    ITstream is(optName, options_[optName]);

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

inline const Foam::string& Foam::argList::operator[](const label index) const
{
    return args_[index];
}


inline const Foam::string& Foam::argList::operator[](const word& optName) const
{
    return options_[optName];
}


// * * * * * * * * * * * * * * * Compatibility * * * * * * * * * * * * * * * //


template<class T>
inline T Foam::argList::argRead(const label index) const
{
    return this->read<T>(index);
}


inline bool Foam::argList::optionFound(const word& optName) const
{
    return found(optName);
}


inline Foam::IStringStream
Foam::argList::optionLookup(const word& optName) const
{
    return lookup(optName);
}


template<class T>
inline T Foam::argList::optionRead(const word& optName) const
{
    return opt<T>(optName);
}


template<class T>
inline bool
Foam::argList::optionReadIfPresent(const word& optName, T& val) const
{
    return readIfPresent<T>(optName, val);
}


template<class T>
inline bool Foam::argList::optionReadIfPresent
(
    const word& optName,
    T& val,
    const T& deflt
) const
{
    return readIfPresent<T>(optName, val, deflt);
}


template<class T>
inline T Foam::argList::optionLookupOrDefault
(
    const word& optName,
    const T& deflt
) const
{
    return lookupOrDefault<T>(optName, deflt);
}


template<class T>
inline Foam::List<T>
Foam::argList::optionReadList(const word& optName) const
{
    return this->readList<T>(optName);
}


// ************************************************************************* //
