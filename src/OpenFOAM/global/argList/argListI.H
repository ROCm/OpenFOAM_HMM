/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2009 OpenCFD Ltd.
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


inline const Foam::fileName& Foam::argList::caseName() const
{
    return case_;
}


inline const Foam::fileName& Foam::argList::globalCaseName() const
{
    return globalCase_;
}


inline Foam::fileName Foam::argList::path() const
{
    return rootPath()/caseName();
}


inline const Foam::stringList& Foam::argList::args() const
{
    return args_;
}


inline const Foam::HashTable<Foam::string>& Foam::argList::options() const
{
    return options_;
}


inline const Foam::string& Foam::argList::option(const word& opt) const
{
    return options_.operator[](opt);
}


inline bool Foam::argList::optionFound(const word& opt) const
{
    return options_.found(opt);
}


inline Foam::IStringStream Foam::argList::optionLookup(const word& opt) const
{
    return IStringStream(option(opt));
}


// * * * * * * * * * * * * Template Specializations  * * * * * * * * * * * * //

namespace Foam
{

    // Template specialization for string
    template<>
    inline Foam::string
    Foam::argList::optionRead<Foam::string>(const word& opt) const
    {
        return option(opt);
    }

    // Template specialization for word
    template<>
    inline Foam::word
    Foam::argList::optionRead<Foam::word>(const word& opt) const
    {
        return option(opt);
    }

    // Template specialization for fileName
    template<>
    inline Foam::fileName
    Foam::argList::optionRead<Foam::fileName>(const word& opt) const
    {
        return option(opt);
    }


} // End namespace Foam


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline T Foam::argList::optionRead(const word& opt) const
{
    T val;

    optionLookup(opt)() >> val;
    return val;
}


template<class T>
inline bool Foam::argList::optionReadIfPresent
(
    const word& opt,
    T& val
) const
{
    if (optionFound(opt))
    {
        val = optionRead<T>(opt);
        return true;
    }
    else
    {
        return false;
    }
}


template<class T>
inline bool Foam::argList::optionReadIfPresent
(
    const word& opt,
    T& val,
    const T& deflt
) const
{
    if (optionReadIfPresent<T>(opt, val))
    {
        return true;
    }
    else
    {
        val = deflt;
        return false;
    }
}


template<class T>
inline T Foam::argList::optionLookupOrDefault
(
    const word& opt,
    const T& deflt
) const
{
    if (optionFound(opt))
    {
        return optionRead<T>(opt);
    }
    else
    {
        return deflt;
    }
}


// ************************************************************************* //
