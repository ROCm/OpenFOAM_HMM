/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2018 OpenCFD Ltd.
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

#include "coordinateSystems.H"
#include "IOPtrList.H"
#include "Time.H"
#include "stringListOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coordinateSystems, 0);
    defineTemplateTypeNameAndDebug(IOPtrList<coordinateSystem>, 0);
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
    // Templated implementation for names() - file-scope
    template<class UnaryMatchPredicate>
    static wordList namesImpl
    (
        const PtrList<coordinateSystem>& list,
        const UnaryMatchPredicate& matcher,
        const bool doSort
    )
    {
        const label len = list.size();

        wordList output(len);

        label count = 0;
        for (label i = 0; i < len; ++i)
        {
            const word& itemName = list[i].name();

            if (matcher(itemName))
            {
                output[count++] = itemName;
            }
        }

        output.resize(count);

        if (doSort)
        {
            Foam::sort(output);
        }

        return output;
    }


    // Templated implementation for indices() - file-scope
    template<class UnaryMatchPredicate>
    static labelList indicesImpl
    (
        const PtrList<coordinateSystem>& list,
        const UnaryMatchPredicate& matcher
    )
    {
        const label len = list.size();

        labelList output(len);

        label count = 0;
        for (label i = 0; i < len; ++i)
        {
            if (matcher(list[i].name()))
            {
                output[count++] = i;
            }
        }

        output.resize(count);

        return output;
    }


    // Templated implementation for findIndex() - file-scope
    template<class UnaryMatchPredicate>
    label findIndexImpl
    (
        const PtrList<coordinateSystem>& list,
        const UnaryMatchPredicate& matcher
    )
    {
        const label len = list.size();

        for (label i = 0; i < len; ++i)
        {
            if (matcher(list[i].name()))
            {
                return i;
            }
        }

        return -1;
    }

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coordinateSystems::coordinateSystems(const IOobject& io)
:
    IOPtrList<coordinateSystem>(io)
{}


Foam::coordinateSystems::coordinateSystems
(
    const IOobject& io,
    const PtrList<coordinateSystem>& content
)
:
    IOPtrList<coordinateSystem>(io, content)
{}


Foam::coordinateSystems::coordinateSystems
(
    const IOobject& io,
    PtrList<coordinateSystem>&& content
)
:
    IOPtrList<coordinateSystem>(io, std::move(content))
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

const Foam::coordinateSystems& Foam::coordinateSystems::New
(
    const objectRegistry& obr
)
{
    // Previously registered?

    const coordinateSystems* ptr =
        obr.lookupObjectPtr<coordinateSystems>(typeName);

    if (ptr)
    {
        return *ptr;
    }

    // Read construct from registry
    return obr.store
    (
        new coordinateSystems
        (
            IOobject
            (
                typeName,
                obr.time().constant(),
                obr,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            )
        )
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::coordinateSystems::indices(const keyType& key) const
{
    if (key.empty())
    {
        return labelList();
    }
    else if (key.isPattern())
    {
        // Match as regex
        regExp matcher(key);
        return indicesImpl(*this, matcher);
    }
    else
    {
        // Compare as literal string
        const word& matcher = key;
        return indicesImpl(*this, matcher);
    }
}


Foam::labelList Foam::coordinateSystems::indices(const wordRes& matcher) const
{
    return indicesImpl(*this, matcher);
}


Foam::label Foam::coordinateSystems::findIndex(const keyType& key) const
{
    if (key.empty())
    {
        return -1;
    }

    if (key.isPattern())
    {
        // Find as regex
        regExp matcher(key);
        return findIndexImpl(*this, matcher);
    }
    else
    {
        // Find as literal string
        const word& matcher = key;
        return findIndexImpl(*this, matcher);
    }
}


Foam::label Foam::coordinateSystems::findIndex(const wordRes& matcher) const
{
    return findIndexImpl(*this, matcher);
}


bool Foam::coordinateSystems::found(const keyType& key) const
{
    return findIndex(key) != -1;
}


Foam::wordList Foam::coordinateSystems::names() const
{
    const PtrList<coordinateSystem>& systems = *this;

    wordList list(systems.size());

    forAll(systems, i)
    {
        list[i] = systems[i].name();
    }

    return list;
}


Foam::wordList Foam::coordinateSystems::names(const keyType& key) const
{
    if (key.empty())
    {
        return wordList();
    }
    else if (key.isPattern())
    {
        // Find as regex
        regExp matcher(key);
        return namesImpl(*this, matcher, false);
    }
    else
    {
        // Find as literal string
        const word& matcher = key;
        return namesImpl(*this, matcher, false);
    }
}


Foam::wordList Foam::coordinateSystems::names(const wordRe& matcher) const
{
    return namesImpl(*this, matcher, false);
}


Foam::wordList Foam::coordinateSystems::names(const wordRes& matcher) const
{
    return namesImpl(*this, matcher, false);
}


bool Foam::coordinateSystems::writeData(Ostream& os) const
{
    os << nl << size() << nl << token::BEGIN_LIST;

    forAll(*this, i)
    {
        os << nl;
        operator[](i).writeDict(os, true);
    }

    os << token::END_LIST << nl;

    return os.good();
}


// ************************************************************************* //
