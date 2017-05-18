/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
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

#include "objectRegistry.H"
#include "stringListOps.H"
#include "predicates.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Templated implementation for classes()
template<class UnaryMatchPredicate>
Foam::HashTable<Foam::wordHashSet> Foam::objectRegistry::classesImpl
(
    const objectRegistry& list,
    const UnaryMatchPredicate& matcher
)
{
    HashTable<wordHashSet> summary(2*list.size());

    // Summary (key,val) = (class-name, object-names)
    forAllConstIters(list, iter)
    {
        if (matcher(iter.key()))
        {
            // Create entry (if needed) and insert
            summary(iter.object()->type()).insert(iter.key());
        }
    }

    return summary;
}


// Templated implementation for names()
template<class Type, class UnaryMatchPredicate>
Foam::wordList Foam::objectRegistry::namesImpl
(
    const objectRegistry& list,
    const UnaryMatchPredicate& matcher,
    const bool doSort
)
{
    wordList objNames(list.size());

    label count = 0;
    forAllConstIters(list, iter)
    {
        if (isA<Type>(*iter()) && matcher(iter()->name()))
        {
            objNames[count++] = iter()->name();
        }
    }

    objNames.setSize(count);

    if (doSort)
    {
        Foam::sort(objNames);
    }

    return objNames;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::wordList Foam::objectRegistry::names() const
{
    return namesImpl<Type>(*this, predicates::always(), false);
}


template<class Type>
Foam::wordList Foam::objectRegistry::names(const wordRe& matcher) const
{
    return namesImpl<Type>(*this, matcher, false);
}


template<class Type>
Foam::wordList Foam::objectRegistry::names
(
    const wordRes& matcher
) const
{
    return namesImpl<Type>(*this, matcher, false);
}


template<class Type>
Foam::wordList Foam::objectRegistry::sortedNames() const
{
    return namesImpl<Type>(*this, predicates::always(), true);
}


template<class Type>
Foam::wordList Foam::objectRegistry::sortedNames
(
    const wordRe& matcher
) const
{
    return namesImpl<Type>(*this, matcher, true);
}


template<class Type>
Foam::wordList Foam::objectRegistry::sortedNames
(
    const wordRes& matcher
) const
{
    return namesImpl<Type>(*this, matcher, true);
}


template<class Type>
Foam::HashTable<const Type*> Foam::objectRegistry::lookupClass
(
    const bool strict
) const
{
    HashTable<const Type*> objectsOfClass(size());

    forAllConstIters(*this, iter)
    {
        if (strict ? isType<Type>(*iter()) : isA<Type>(*iter()))
        {
            objectsOfClass.insert
            (
                iter()->name(),
                dynamic_cast<const Type*>(iter())
            );
        }
    }

    return objectsOfClass;
}


template<class Type>
Foam::HashTable<Type*> Foam::objectRegistry::lookupClass
(
    const bool strict
)
{
    HashTable<Type*> objectsOfClass(size());

    forAllIters(*this, iter)
    {
        if (strict ? isType<Type>(*iter()) : isA<Type>(*iter()))
        {
            objectsOfClass.insert
            (
                iter()->name(),
                dynamic_cast<Type*>(iter())
            );
        }
    }

    return objectsOfClass;
}


template<class Type>
bool Foam::objectRegistry::foundObject
(
    const word& name,
    const bool recursive
) const
{
    const Type* ptr = this->lookupObjectPtr<Type>(name, recursive);

    if (ptr)
    {
        return true;
    }
    else
    {
        return false;
    }
}


template<class Type>
const Type& Foam::objectRegistry::lookupObject
(
    const word& name,
    const bool recursive
) const
{
    const_iterator iter = find(name);

    if (iter.found())
    {
        const Type* ptr = dynamic_cast<const Type*>(iter());

        if (ptr)
        {
            return *ptr;
        }

        FatalErrorInFunction
            << nl
            << "    lookup of " << name << " from objectRegistry "
            << this->name()
            << " successful\n    but it is not a " << Type::typeName
            << ", it is a " << iter()->type()
            << abort(FatalError);
    }
    else if (recursive && this->parentNotTime())
    {
        return parent_.lookupObject<Type>(name, recursive);
    }

    FatalErrorInFunction
        << nl
        << "    request for " << Type::typeName
        << " " << name << " from objectRegistry " << this->name()
        << " failed\n    available objects of type " << Type::typeName
        << " are" << nl
        << names<Type>()
        << abort(FatalError);

    return NullObjectRef<Type>();
}


template<class Type>
Type& Foam::objectRegistry::lookupObjectRef
(
    const word& name,
    const bool recursive
) const
{
    const Type& ref = this->lookupObject<Type>(name, recursive);
    // The above will already fail if things didn't work

    return const_cast<Type&>(ref);
}


template<class Type>
const Type* Foam::objectRegistry::lookupObjectPtr
(
    const word& name,
    const bool recursive
) const
{
    const_iterator iter = find(name);

    if (iter.found())
    {
        const Type* ptr = dynamic_cast<const Type*>(iter());

        if (ptr)
        {
            return ptr;
        }
    }
    else if (recursive && this->parentNotTime())
    {
        return parent_.lookupObjectPtr<Type>(name, recursive);
    }

    return nullptr;
}


template<class Type>
Type* Foam::objectRegistry::lookupObjectRefPtr
(
    const word& name,
    const bool recursive
) const
{
    const Type* ptr = this->lookupObjectPtr<Type>(name, recursive);

    return const_cast<Type*>(ptr);
}


// ************************************************************************* //
