/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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

Namespace
    Foam::predicates

Description
    Various constant predicate types.

SourceFiles
    predicates.H

\*---------------------------------------------------------------------------*/

#ifndef predicates_H
#define predicates_H

#include <string>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace predicates
{

/*---------------------------------------------------------------------------*\
                           Class always Declaration
\*---------------------------------------------------------------------------*/

//- Unary and binary predicates that always return true, useful for templating.
struct always
{
    typedef always value_type;

    //- Null constructible
    inline always()
    {}

    //- Evaluated as a bool
    //  \return true
    inline operator bool() const
    {
        return true;
    }

    //- Unary predicate
    //  \return true
    template<class T>
    inline bool operator()(const T&) const
    {
        return true;
    }

    //- Binary predicate
    //  \return true
    template<class T1, class T2>
    inline bool operator()(const T1&, const T2&) const
    {
        return true;
    }

    //- String match
    //  \return true
    inline bool match(const std::string&, bool literal=false) const
    {
        return true;
    }
};


/*---------------------------------------------------------------------------*\
                            Class never Declaration
\*---------------------------------------------------------------------------*/

//- Unary and binary predicates that never return true, useful for templating.
struct never
{
    typedef never value_type;

    //- Null constructible
    inline never()
    {}

    //- Evaluated as a bool
    //  \return false
    inline operator bool() const
    {
        return false;
    }

    //- Unary predicate
    //  \return false
    template<class T>
    inline bool operator()(const T&) const
    {
        return false;
    }

    //- Binary predicate
    //  \return false
    template<class T1, class T2>
    inline bool operator()(const T1&, const T2&) const
    {
        return false;
    }

    //- String match
    //  \return false
    inline bool match(const std::string&, bool literal=false) const
    {
        return false;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace predicates

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
