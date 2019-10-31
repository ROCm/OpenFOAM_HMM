/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
    Copyright (C) 2018-2019 OpenCFD Ltd.
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

Class
    Foam::Hash

Description
    Hash function class.
    The default definition is for primitives, non-primitives used to hash
    entries on hash tables likely need a specialized version of this class.

\*---------------------------------------------------------------------------*/

#ifndef Hash_H
#define Hash_H

#include "label.H"
#include "uLabel.H"
#include "Hasher.H"
#include "fileName.H"
#include "wordRe.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class Hash Declaration
\*---------------------------------------------------------------------------*/

template<class T>
struct Hash
{
    inline unsigned operator()(const T& obj, unsigned seed=0) const
    {
        return Hasher(&obj, sizeof(obj), seed);
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Hash specialization for label
template<>
struct Hash<Foam::label>
{
    //- Incrementally hash a label.
    //  This will necessarily return a different value than the
    //  non-incremental version.
    inline unsigned operator()(const label obj, unsigned seed) const
    {
        return Hasher(&obj, sizeof(label), seed);
    }

    //- Return the unsigned representation of a label.
    //  This helps if people have relied on the hash value corresponding to
    //  the natural order.
    inline unsigned operator()(const label obj) const
    {
        return obj;
    }
};


//- Hash specialization for string
template<>
struct Hash<Foam::string>
{
    inline unsigned operator()(const string& obj, unsigned seed=0) const
    {
        return string::hash()(obj, seed);
    }
};


//- Hash specialization for word
template<>
struct Hash<Foam::word>
{
    inline unsigned operator()(const word& obj, unsigned seed=0) const
    {
        return string::hash()(obj, seed);
    }
};


//- Hash specialization for fileName
template<>
struct Hash<Foam::fileName>
{
    inline unsigned operator()(const fileName& obj, unsigned seed=0) const
    {
        return string::hash()(obj, seed);
    }
};


//- Hash specialization for wordRe
template<>
struct Hash<Foam::wordRe>
{
    inline unsigned operator()(const wordRe& obj, unsigned seed=0) const
    {
        return string::hash()(obj, seed);
    }
};


//- Hash specialization for keyType
template<>
struct Hash<Foam::keyType>
{
    inline unsigned operator()(const keyType& obj, unsigned seed=0) const
    {
        return string::hash()(obj, seed);
    }
};


//- Hash specialization for pointers, interpret pointer as a integer type.
template<>
struct Hash<void*>
{
    inline unsigned operator()(const void* const& obj, unsigned seed) const
    {
        return Hash<intptr_t>()(intptr_t(obj), seed);
    }

    inline unsigned operator()(const void* const& obj) const
    {
        return Hash<intptr_t>()(intptr_t(obj));
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
