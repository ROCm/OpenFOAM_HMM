/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeName(coordinateSystems);
}

// File-local

//- Header name for 1806 and earlier
static const char* headerTypeCompat = "IOPtrList<coordinateSystem>";


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


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::coordinateSystems::readFromStream(const bool valid)
{
    Istream& is = readStream(word::null, valid);

    if (valid)
    {
        if (headerClassName() == typeName)
        {
            this->readIstream(is, coordinateSystem::iNew());
            close();
        }
        else if (headerClassName() == headerTypeCompat)
        {
            // Older (1806 and earlier) header name
            std::cerr
                << "--> FOAM IOWarning :" << nl
                << "    Found header class name '" << headerTypeCompat
                << "' instead of '" << typeName << "'" << nl;

            error::warnAboutAge("header class", 1806);

            this->readIstream(is, coordinateSystem::iNew());
            close();
        }
        else
        {
            FatalIOErrorInFunction(is)
                << "unexpected class name " << headerClassName()
                << " expected " << typeName
                << " or " << headerTypeCompat << nl
                << "    while reading object " << name()
                << exit(FatalIOError);
        }
    }
}


bool Foam::coordinateSystems::readObject(const IOobject& io)
{
    if
    (
        (
            io.readOpt() == IOobject::MUST_READ
         || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
        )
     || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        readFromStream();
        return true;
    }

    return false;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coordinateSystems::coordinateSystems(const IOobject& io)
:
    regIOobject(io),
    PtrList<coordinateSystem>()
{
    readObject(io);
}


Foam::coordinateSystems::coordinateSystems
(
    const IOobject& io,
    const PtrList<coordinateSystem>& content
)
:
    regIOobject(io),
    PtrList<coordinateSystem>()
{
    if (!readObject(io))
    {
        static_cast<PtrList<coordinateSystem>&>(*this) = content;
    }
}


Foam::coordinateSystems::coordinateSystems
(
    const IOobject& io,
    PtrList<coordinateSystem>&& content
)
:
    regIOobject(io),
    PtrList<coordinateSystem>(std::move(content))
{
    readObject(io);
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

const Foam::coordinateSystems& Foam::coordinateSystems::New
(
    const objectRegistry& obr
)
{
    // Previously registered?

    const coordinateSystems* ptr = obr.findObject<coordinateSystems>(typeName);

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


const Foam::coordinateSystem*
Foam::coordinateSystems::lookupPtr(const word& name) const
{
    const label index = this->findIndex(name);

    if (coordinateSystem::debug)
    {
        InfoInFunction
            << "Global coordinate system: "
            << name << "=" << index << endl;
    }

    if (index < 0)
    {
        return nullptr;
    }

    return this->operator()(index);
}


const Foam::coordinateSystem&
Foam::coordinateSystems::lookup(const word& name) const
{
    const label index = this->findIndex(name);

    if (index < 0)
    {
        FatalErrorInFunction
            << "Could not find coordinate system: " << name << nl
            << "available coordinate systems: "
            << flatOutput(names()) << nl << nl
            << exit(FatalError);
    }
    if (coordinateSystem::debug)
    {
        InfoInFunction
            << "Global coordinate system: "
            << name << "=" << index << endl;
    }

    return this->operator[](index);
}


Foam::wordList Foam::coordinateSystems::names() const
{
    const PtrList<coordinateSystem>& list = *this;

    wordList result(list.size());

    forAll(list, i)
    {
        result[i] = list[i].name();
    }

    return result;
    // return ListOps::create<word>(list, nameOp<coordinateSystem>());
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
    const PtrList<coordinateSystem>& list = *this;

    os << nl << size() << nl << token::BEGIN_LIST;

    for (const coordinateSystem& csys : list)
    {
        os << nl;
        csys.writeEntry(csys.name(), os);
    }

    os << token::END_LIST << nl;

    return os.good();
}


bool Foam::coordinateSystems::writeObject
(
    IOstream::streamFormat,
    IOstream::versionNumber ver,
    IOstream::compressionType,
    const bool valid
) const
{
    // Force ASCII writing
    return regIOobject::writeObject
    (
        IOstream::ASCII,
        ver,
        IOstream::UNCOMPRESSED,
        valid
    );
}


// ************************************************************************* //
