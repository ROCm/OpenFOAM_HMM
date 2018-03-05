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
        const IOPtrList<coordinateSystem>& list,
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

Foam::labelList Foam::coordinateSystems::findIndices(const keyType& key) const
{
    labelList indices;
    if (key.empty())
    {
        // no-op
    }
    else if (key.isPattern())
    {
        indices = findStrings(key, this->names());
    }
    else
    {
        indices.setSize(this->size());
        label count = 0;
        forAll(*this, i)
        {
            if (key == operator[](i).name())
            {
                indices[count++] = i;
            }
        }
        indices.setSize(count);
    }

    return indices;
}


Foam::label Foam::coordinateSystems::findIndex(const keyType& key) const
{
    if (key.empty())
    {
        // no-op
    }
    else if (key.isPattern())
    {
        labelList indices = findIndices(key);
        if (!indices.empty())
        {
            return indices.first();  // first match
        }
    }
    else
    {
        forAll(*this, i)
        {
            if (key == operator[](i).name())
            {
                return i;
            }
        }
    }

    // Not found
    return -1;
}


bool Foam::coordinateSystems::found(const keyType& key) const
{
    return findIndex(key) != -1;
}


Foam::wordList Foam::coordinateSystems::names() const
{
    wordList output(size());

    forAll(*this, i)
    {
        output[i] = operator[](i).name();
    }

    return output;
}


Foam::wordList Foam::coordinateSystems::names(const keyType& matcher) const
{
    return
    (
        matcher.isPattern()
      ? namesImpl(*this, regExp(matcher), false)
      : namesImpl(*this, matcher, false)
    );
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
