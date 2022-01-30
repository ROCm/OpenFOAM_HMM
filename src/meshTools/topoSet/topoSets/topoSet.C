/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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

#include "topoSet.H"
#include "mapPolyMesh.H"
#include "polyMesh.H"
#include "boundBox.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(topoSet, 0);
    defineRunTimeSelectionTable(topoSet, word);
    defineRunTimeSelectionTable(topoSet, size);
    defineRunTimeSelectionTable(topoSet, set);

    int Foam::topoSet::disallowGenericSets
    (
        debug::debugSwitch("disallowGenericSets", 0)
    );
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::topoSet>
Foam::topoSet::New
(
    const word& setType,
    const polyMesh& mesh,
    const word& name,
    readOption r,
    writeOption w
)
{
    auto* ctorPtr = wordConstructorTable(setType);

    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            "set",
            setType,
            *wordConstructorTablePtr_
        ) << exit(FatalError);
    }

    return autoPtr<topoSet>(ctorPtr(mesh, name, r, w));
}


Foam::autoPtr<Foam::topoSet>
Foam::topoSet::New
(
    const word& setType,
    const polyMesh& mesh,
    const word& name,
    const label size,
    writeOption w
)
{
    auto* ctorPtr = sizeConstructorTable(setType);

    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            "set",
            setType,
            *sizeConstructorTablePtr_
        ) << exit(FatalError);
    }

    return autoPtr<topoSet>(ctorPtr(mesh, name, size, w));
}


Foam::autoPtr<Foam::topoSet>
Foam::topoSet::New
(
    const word& setType,
    const polyMesh& mesh,
    const word& name,
    const topoSet& set,
    writeOption w
)
{
    auto* ctorPtr = setConstructorTable(setType);

    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            "set",
            setType,
            *setConstructorTablePtr_
        ) << exit(FatalError);
    }

    return autoPtr<topoSet>(ctorPtr(mesh, name, set, w));
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::fileName Foam::topoSet::localPath
(
    const polyMesh& mesh,
    const word& name
)
{
    return mesh.facesInstance()/mesh.dbDir()/polyMesh::meshSubDir/"sets"/name;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

// Update stored cell numbers using map.
// Do in two passes to prevent allocation if nothing changed.
void Foam::topoSet::updateLabels(const labelUList& map)
{
    labelHashSet& labels = *this;

    // Iterate over map to see if anything changed

    bool changed = false;

    for (const label oldId : labels)
    {
        if (oldId < 0 || oldId >= map.size())
        {
            FatalErrorInFunction
                << "Illegal content " << oldId << " of set:" << name()
                << " of type " << type() << nl
                << "Value should be between [0," << map.size() << ')'
                << endl
                << abort(FatalError);
        }

        const label newId = map[oldId];

        if (newId != oldId)
        {
            changed = true;
            #ifdef FULLDEBUG
            continue;  // Check all elements in FULLDEBUG mode
            #endif
            break;
        }
    }

    if (!changed)
    {
        return;
    }


    // Relabel. Use second labelHashSet to prevent overlapping.

    labelHashSet newLabels(2*labels.size());

    for (const label oldId : labels)
    {
        const label newId = map[oldId];

        if (newId >= 0)
        {
            newLabels.set(newId);
        }
    }

    labels.transfer(newLabels);
}


void Foam::topoSet::check(const label maxSize)
{
    const labelHashSet& labels = *this;

    for (const label oldId : labels)
    {
        if (oldId < 0 || oldId >= maxSize)
        {
            FatalErrorInFunction
                << "Illegal content " << oldId << " of set:" << name()
                << " of type " << type() << nl
                << "Value should be between [0," << maxSize << ')'
                << endl
                << abort(FatalError);
        }
    }
}


// Write maxElem elements, starting at iter. Updates iter and elemI.
void Foam::topoSet::writeDebug
(
    Ostream& os,
    const label maxElem,
    topoSet::const_iterator& iter,
    label& elemI
) const
{
    label n = 0;

    for (; (iter != cend()) && (n < maxElem); ++iter)
    {
        if (n && ((n % 10) == 0))
        {
            os << nl;
        }
        os << iter.key() << ' ';

        ++n;
        ++elemI;
    }
}


// Write maxElem elements, starting at iter. Updates iter and elemI.
void Foam::topoSet::writeDebug
(
    Ostream& os,
    const pointField& coords,
    const label maxElem,
    topoSet::const_iterator& iter,
    label& elemI
) const
{
    label n = 0;

    for (; (iter != cend()) && (n < maxElem); ++iter)
    {
        if (n && ((n % 3) == 0))
        {
            os << nl;
        }
        os << iter.key() << coords[iter.key()] << ' ';

        ++n;
        ++elemI;
    }
}


void Foam::topoSet::writeDebug
(
    Ostream& os,
    const pointField& coords,
    const label maxLen
) const
{
    // Bounding box of contents.
    boundBox bb(pointField(coords, toc()), true);

    os  << "Set bounding box: min = "
        << bb.min() << "    max = " << bb.max() << " metres." << nl << endl;

    label n = 0;

    topoSet::const_iterator iter = this->cbegin();

    if (size() <= maxLen)
    {
        writeDebug(os, coords, maxLen, iter, n);
    }
    else
    {
        label halfLen = maxLen/2;

        os  << "Size larger than " << maxLen << ". Printing first and last "
            << halfLen << " elements:" << nl << endl;

        writeDebug(os, coords, halfLen, iter, n);

        os  << nl << "  .." << nl << endl;

        for (; n < size() - halfLen; ++n)
        {
            ++iter;
        }

        writeDebug(os, coords, halfLen, iter, n);
    }
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::IOobject Foam::topoSet::findIOobject
(
    const polyMesh& mesh,
    const word& name,
    readOption r,
    writeOption w
)
{
    IOobject io
    (
        name,
        mesh.time().findInstance
        (
            mesh.dbDir()/polyMesh::meshSubDir/"sets",
            word::null,
            IOobject::READ_IF_PRESENT,
            mesh.facesInstance()
        ),
        polyMesh::meshSubDir/"sets",
        mesh,
        r,
        w
    );

    if (!io.typeHeaderOk<topoSet>(false) && disallowGenericSets != 0)
    {
        DebugInfo<< "Setting no read for set " << name << endl;
        io.readOpt(IOobject::NO_READ);
    }

    return io;
}


Foam::IOobject Foam::topoSet::findIOobject
(
    const Time& runTime,
    const word& name,
    readOption r,
    writeOption w
)
{
    return IOobject
    (
        name,
        runTime.findInstance
        (
            polyMesh::meshSubDir/"sets",
            word::null,
            IOobject::MUST_READ,
            runTime.findInstance
            (
                polyMesh::meshSubDir,
                "faces",
                IOobject::READ_IF_PRESENT
            )
        ),
        polyMesh::meshSubDir/"sets",
        runTime,
        r,
        w
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::topoSet::topoSet(const IOobject& obj, const word& wantedType)
:
    regIOobject(obj)
{
    if
    (
        readOpt() == IOobject::MUST_READ
     || readOpt() == IOobject::MUST_READ_IF_MODIFIED
     || (readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        if (readStream(wantedType).good())
        {
            readStream(wantedType) >> static_cast<labelHashSet&>(*this);

            close();
        }
    }
}


Foam::topoSet::topoSet
(
    const polyMesh& mesh,
    const word& wantedType,
    const word& name,
    readOption r,
    writeOption w
)
:
    regIOobject(findIOobject(mesh, name, r, w))
{
    if
    (
        readOpt() == IOobject::MUST_READ
     || readOpt() == IOobject::MUST_READ_IF_MODIFIED
     || (readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        if (readStream(wantedType).good())
        {
            readStream(wantedType) >> static_cast<labelHashSet&>(*this);

            close();
        }
    }
}


Foam::topoSet::topoSet
(
    const polyMesh& mesh,
    const word& name,
    const label size,
    writeOption w
)
:
    regIOobject(findIOobject(mesh, name, IOobject::NO_READ, w)),
    labelHashSet(size)
{}


Foam::topoSet::topoSet
(
    const polyMesh& mesh,
    const word& name,
    const labelHashSet& labels,
    writeOption w
)
:
    regIOobject(findIOobject(mesh, name, IOobject::NO_READ, w)),
    labelHashSet(labels)
{}


Foam::topoSet::topoSet
(
    const polyMesh& mesh,
    const word& name,
    labelHashSet&& labels,
    writeOption w
)
:
    regIOobject(findIOobject(mesh, name, IOobject::NO_READ, w)),
    labelHashSet(std::move(labels))
{}


Foam::topoSet::topoSet
(
    const polyMesh& mesh,
    const word& name,
    const labelUList& labels,
    writeOption w
)
:
    regIOobject(findIOobject(mesh, name, IOobject::NO_READ, w)),
    labelHashSet(labels)
{}


Foam::topoSet::topoSet(const IOobject& io, const label size)
:
    regIOobject(io),
    labelHashSet(size)
{}


Foam::topoSet::topoSet(const IOobject& io, const labelHashSet& labels)
:
    regIOobject(io),
    labelHashSet(labels)
{}


Foam::topoSet::topoSet(const IOobject& io, labelHashSet&& labels)
:
    regIOobject(io),
    labelHashSet(std::move(labels))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::topoSet::found(const label id) const
{
    return static_cast<const labelHashSet&>(*this).found(id);
}


bool Foam::topoSet::set(const label id)
{
    return static_cast<labelHashSet&>(*this).set(id);
}


bool Foam::topoSet::unset(const label id)
{
    return static_cast<labelHashSet&>(*this).unset(id);
}


void Foam::topoSet::set(const labelUList& labels)
{
    static_cast<labelHashSet&>(*this).set(labels);
}


void Foam::topoSet::unset(const labelUList& labels)
{
    static_cast<labelHashSet&>(*this).unset(labels);
}


void Foam::topoSet::invert(const label maxLen)
{
    // Retain a copy of the original (current) set.
    labelHashSet original
    (
        std::move(static_cast<labelHashSet&>(*this))
    );

    clear();  // Maybe don't trust the previous move operation
    resize(2*max(64, (maxLen - original.size())));

    for (label id=0; id < maxLen; ++id)
    {
        if (!original.found(id))
        {
            this->set(id);
        }
    }
}


void Foam::topoSet::subset(const topoSet& set)
{
    // Only retain entries found in both sets
    static_cast<labelHashSet&>(*this) &= set;
}


void Foam::topoSet::addSet(const topoSet& set)
{
    // Add entries to the set
    static_cast<labelHashSet&>(*this) |= set;
}


void Foam::topoSet::subtractSet(const topoSet& set)
{
    // Subtract entries from the set
    static_cast<labelHashSet&>(*this) -= set;
}


void Foam::topoSet::deleteSet(const topoSet& set)
{
    this->subtractSet(set);
}


void Foam::topoSet::sync(const polyMesh&)
{
    NotImplemented;
}


void Foam::topoSet::writeDebug(Ostream& os, const label maxLen) const
{
    label n = 0;

    topoSet::const_iterator iter = this->cbegin();

    if (size() <= maxLen)
    {
        writeDebug(os, maxLen, iter, n);
    }
    else
    {
        label halfLen = maxLen/2;

        os  << "Size larger than " << maxLen << ". Printing first and last "
            << halfLen << " elements:" << nl << endl;

        writeDebug(os, halfLen, iter, n);

        os  << nl << "  .." << nl << endl;

        for (; n < size() - halfLen; ++n)
        {
            ++iter;
        }

        writeDebug(os, halfLen, iter, n);
    }
}


bool Foam::topoSet::writeData(Ostream& os) const
{
    return (os << *this).good();
}


void Foam::topoSet::updateMesh(const mapPolyMesh&)
{
    NotImplemented;
}


void Foam::topoSet::removeFiles(const polyMesh& mesh)
{
    IOobject io
    (
        "dummy",
        mesh.facesInstance(),
        mesh.meshSubDir/"sets",
        mesh
    );
    fileName setsDir(io.path());

    if (debug) DebugVar(setsDir);

    if (isDir(setsDir))
    {
        rmDir(setsDir);
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::topoSet::operator=(const topoSet& rhs)
{
    labelHashSet::operator=(rhs);
}


// ************************************************************************* //
