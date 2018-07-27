/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2018 OpenCFD Ltd.
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

#include "IOobjectList.H"
#include "Time.H"
#include "OSspecific.H"
#include "IOList.H"
#include "predicates.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{
    // Templated implementation for lookup() - file-scope
    template<class UnaryMatchPredicate>
    static IOobjectList lookupImpl
    (
        const IOobjectList& list,
        const UnaryMatchPredicate& matcher
    )
    {
        IOobjectList results(list.size());

        forAllConstIters(list, iter)
        {
            const word& key = iter.key();
            const IOobject* io = iter.object();

            if (matcher(key))
            {
                if (IOobject::debug)
                {
                    InfoInFunction << "Found " << key << endl;
                }

                results.set(key, new IOobject(*io));
            }
        }

        return results;
    }


    // Templated implementation for lookupClass() - file-scope
    template<class UnaryMatchPredicate>
    static IOobjectList lookupClassImpl
    (
        const IOobjectList& list,
        const word& clsName,
        const UnaryMatchPredicate& matcher
    )
    {
        IOobjectList results(list.size());

        forAllConstIters(list, iter)
        {
            const word& key = iter.key();
            const IOobject* io = iter.object();

            if (clsName == io->headerClassName() && matcher(key))
            {
                if (IOobject::debug)
                {
                    InfoInFunction << "Found " << key << endl;
                }

                results.set(key, new IOobject(*io));
            }
        }

        return results;
    }


    // Templated implementation for classes() - file-scope
    template<class UnaryMatchPredicate>
    static HashTable<wordHashSet> classesImpl
    (
        const IOobjectList& list,
        const UnaryMatchPredicate& matcher
    )
    {
        HashTable<wordHashSet> summary(2*list.size());

        // Summary (key,val) = (class-name, object-names)
        forAllConstIters(list, iter)
        {
            const word& key = iter.key();
            const IOobject* io = iter.object();

            if (matcher(key))
            {
                // Create entry (if needed) and insert
                summary(io->headerClassName()).insert(key);
            }
        }

        return summary;
    }


    // Templated implementation for names(), sortedNames() - file-scope
    template<class UnaryMatchPredicate>
    static wordList namesImpl
    (
        const IOobjectList& list,
        const word& clsName,
        const UnaryMatchPredicate& matcher,
        const bool doSort
    )
    {
        wordList objNames(list.size());

        label count = 0;
        forAllConstIters(list, iter)
        {
            const word& key = iter.key();
            const IOobject* io = iter.object();

            if (clsName == io->headerClassName() && matcher(key))
            {
                objNames[count] = key;
                ++count;
            }
        }

        objNames.resize(count);

        if (doSort)
        {
            Foam::sort(objNames);
        }

        return objNames;
    }


    // With syncPar = true, check that object names are the same on
    // all processors. Trigger FatalError if not.
    //
    // The object names are sorted as a side-effect, since this is
    // required for consistent ordering across all processors.
    static bool checkNames(wordList& masterNames, const bool syncPar)
    {
        Foam::sort(masterNames);

        if (syncPar && Pstream::parRun())
        {
            const wordList localNames(masterNames);
            Pstream::scatter(masterNames);

            if (localNames != masterNames)
            {
                FatalErrorInFunction
                    << "Objects not synchronised across processors." << nl
                    << "Master has " << flatOutput(masterNames) << nl
                    << "Processor " << Pstream::myProcNo()
                    << " has " << flatOutput(localNames)
                    << exit(FatalError);

                return false;
            }
        }

        return true;
    }

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IOobjectList::IOobjectList()
:
    HashPtrTable<IOobject>()
{}


Foam::IOobjectList::IOobjectList(const label nObjects)
:
    HashPtrTable<IOobject>(nObjects)  // Could also use 2*nObjects instead
{}


Foam::IOobjectList::IOobjectList(const IOobjectList& list)
:
    HashPtrTable<IOobject>(list)
{}


Foam::IOobjectList::IOobjectList(IOobjectList&& list)
:
    HashPtrTable<IOobject>(std::move(list))
{}


Foam::IOobjectList::IOobjectList
(
    const objectRegistry& db,
    const fileName& instance,
    const fileName& local,
    IOobject::readOption r,
    IOobject::writeOption w,
    bool registerObject
)
:
    HashPtrTable<IOobject>()
{
    word newInstance;
    fileNameList objNames = fileHandler().readObjects
    (
        db,
        instance,
        local,
        newInstance
    );

    for (const auto& objName : objNames)
    {
        auto objectPtr = autoPtr<IOobject>::New
        (
            objName,
            newInstance,
            local,
            db,
            r,
            w,
            registerObject
        );

        bool ok = false;
        const bool throwingIOerr = FatalIOError.throwExceptions();

        try
        {
            // Use object with local scope and current instance (no searching)
            ok = objectPtr->typeHeaderOk<IOList<label>>(false, false);
        }
        catch (Foam::IOerror& err)
        {
            Warning
                << err << nl << endl;
        }

        FatalIOError.throwExceptions(throwingIOerr);

        if (ok)
        {
            insert(objectPtr->name(), objectPtr);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::IOobjectList::add(autoPtr<IOobject>& objectPtr)
{
    if (objectPtr.valid())
    {
        return insert(objectPtr->name(), objectPtr);
    }

    return false;
}


bool Foam::IOobjectList::add(autoPtr<IOobject>&& objectPtr)
{
    if (objectPtr.valid())
    {
        return insert(objectPtr->name(), objectPtr);
    }

    return false;
}


Foam::label Foam::IOobjectList::append(const IOobjectList& other)
{
    label count = 0;

    forAllConstIters(other, iter)
    {
        if (!found(iter.key()))
        {
            if (IOobject::debug)
            {
                InfoInFunction << "Copy append " << iter.key() << nl;
            }

            set(iter.key(), new IOobject(*(iter.object())));
            ++count;
        }
    }

    return count;
}


Foam::label Foam::IOobjectList::append(IOobjectList&& other)
{
    // Remove by name to avoid uncertainties about invalid iterators

    label count = 0;

    wordList keys(other.toc());

    for (const word& key : keys)
    {
        if (!found(key))
        {
            if (IOobject::debug)
            {
                InfoInFunction << "Move append " << key << nl;
            }

            if (add(other.remove(key)))
            {
                ++count;
            }
        }
    }

    return count;
}


bool Foam::IOobjectList::remove(const IOobject& io)
{
    return erase(io.name());
}


Foam::IOobject* Foam::IOobjectList::lookup(const word& name) const
{
    const_iterator iter = cfind(name);

    if (iter.found())
    {
        if (IOobject::debug)
        {
            InfoInFunction << "Found " << name << endl;
        }

        return const_cast<IOobject*>(*iter);
    }

    if (IOobject::debug)
    {
        InfoInFunction << "Could not find " << name << endl;
    }
    return nullptr;
}


Foam::IOobjectList Foam::IOobjectList::lookup(const wordRe& matcher) const
{
    return lookupImpl(*this, matcher);
}


Foam::IOobjectList Foam::IOobjectList::lookup(const wordRes& matcher) const
{
    return lookupImpl(*this, matcher);
}


Foam::IOobjectList Foam::IOobjectList::lookup(const wordHashSet& matcher) const
{
    return lookupImpl(*this, matcher);
}


Foam::IOobjectList Foam::IOobjectList::lookupClass(const word& clsName) const
{
    return lookupClassImpl(*this, clsName, predicates::always());
}


Foam::HashTable<Foam::wordHashSet> Foam::IOobjectList::classes() const
{
    return classesImpl(*this, predicates::always());
}


Foam::HashTable<Foam::wordHashSet>
Foam::IOobjectList::classes(const wordRe& matcher) const
{
    return classesImpl(*this, matcher);
}


Foam::HashTable<Foam::wordHashSet>
Foam::IOobjectList::classes(const wordRes& matcher) const
{
    return classesImpl(*this, matcher);
}


Foam::HashTable<Foam::wordHashSet>
Foam::IOobjectList::classes(const wordHashSet& matcher) const
{
    return classesImpl(*this, matcher);
}


Foam::wordList Foam::IOobjectList::names() const
{
    return HashPtrTable<IOobject>::toc();
}


Foam::wordList Foam::IOobjectList::sortedNames() const
{
    return HashPtrTable<IOobject>::sortedToc();
}


Foam::wordList Foam::IOobjectList::names(const bool syncPar) const
{
    wordList objNames(HashPtrTable<IOobject>::toc());

    checkNames(objNames, syncPar);
    return objNames;
}


Foam::wordList Foam::IOobjectList::names
(
    const word& clsName
) const
{
    return namesImpl(*this, clsName, predicates::always(), false);
}


Foam::wordList Foam::IOobjectList::sortedNames
(
    const word& clsName
) const
{
    return namesImpl(*this, clsName, predicates::always(), true);
}


Foam::wordList Foam::IOobjectList::names
(
    const word& clsName,
    const bool syncPar
) const
{
    wordList objNames(namesImpl(*this, clsName, predicates::always(), false));

    checkNames(objNames, syncPar);
    return objNames;
}


Foam::wordList Foam::IOobjectList::names
(
    const word& clsName,
    const wordRe& matcher
) const
{
    return namesImpl(*this, clsName, matcher, false);
}


Foam::wordList Foam::IOobjectList::sortedNames
(
    const word& clsName,
    const wordRe& matcher
) const
{
    return namesImpl(*this, clsName, matcher, true);
}


Foam::wordList Foam::IOobjectList::names
(
    const word& clsName,
    const wordRe& matcher,
    const bool syncPar
) const
{
    wordList objNames(namesImpl(*this, clsName, matcher, false));

    checkNames(objNames, syncPar);
    return objNames;
}


Foam::wordList Foam::IOobjectList::names
(
    const word& clsName,
    const wordRes& matcher
) const
{
    return namesImpl(*this, clsName, matcher, false);
}


Foam::wordList Foam::IOobjectList::sortedNames
(
    const word& clsName,
    const wordRes& matcher
) const
{
    return namesImpl(*this, clsName, matcher, true);
}


Foam::wordList Foam::IOobjectList::names
(
    const word& clsName,
    const wordRes& matcher,
    const bool syncPar
) const
{
    wordList objNames(namesImpl(*this, clsName, matcher, false));

    checkNames(objNames, syncPar);
    return objNames;
}


Foam::wordList Foam::IOobjectList::names
(
    const word& clsName,
    const wordHashSet& matcher
) const
{
    return namesImpl(*this, clsName, matcher, false);
}


Foam::wordList Foam::IOobjectList::sortedNames
(
    const word& clsName,
    const wordHashSet& matcher
) const
{
    return namesImpl(*this, clsName, matcher, true);
}


Foam::wordList Foam::IOobjectList::names
(
    const word& clsName,
    const wordHashSet& matcher,
    const bool syncPar
) const
{
    wordList objNames(namesImpl(*this, clsName, matcher, false));

    checkNames(objNames, syncPar);
    return objNames;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::IOobjectList::operator=(IOobjectList&& list)
{
    transfer(list);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const IOobjectList& list)
{
    os << nl << list.size() << nl << token::BEGIN_LIST << nl;

    forAllConstIters(list, it)
    {
        os << it.key() << token::SPACE << it.object()->headerClassName() << nl;
    }

    os << token::END_LIST;
    os.check(FUNCTION_NAME);

    return os;
}


// ************************************************************************* //
