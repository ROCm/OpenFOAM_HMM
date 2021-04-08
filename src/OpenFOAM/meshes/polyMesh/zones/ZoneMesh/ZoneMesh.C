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

#include "ZoneMesh.H"
#include "entry.H"
#include "demandDrivenData.H"
#include "Pstream.H"
#include "PtrListOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    template<class ZoneType, class MeshType>
    int Foam::ZoneMesh<ZoneType, MeshType>::disallowGenericZones
    (
        debug::debugSwitch("disallowGenericZones", 0)
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ZoneType, class MeshType>
void Foam::ZoneMesh<ZoneType, MeshType>::calcZoneMap() const
{
    if (zoneMapPtr_)
    {
        FatalErrorInFunction
            << "zone map already calculated"
            << abort(FatalError);
    }
    else
    {
        // Count number of objects in all zones
        label nObjects = 0;

        const PtrList<ZoneType>& zones = *this;

        for (const ZoneType& zn : zones)
        {
            nObjects += zn.size();
        }

        zoneMapPtr_ = new Map<label>(2*nObjects);
        Map<label>& zm = *zoneMapPtr_;

        // Fill in objects of all zones into the map.
        // The key is the global object index, value is the zone index

        label zonei = 0;

        for (const ZoneType& zn : zones)
        {
            const labelList& labels = zn;

            for (const label idx : labels)
            {
                zm.insert(idx, zonei);
            }

            ++zonei;
        }
    }
}


template<class ZoneType, class MeshType>
bool Foam::ZoneMesh<ZoneType, MeshType>::read()
{
    if
    (
        readOpt() == IOobject::MUST_READ
     || readOpt() == IOobject::MUST_READ_IF_MODIFIED
     || (readOpt() == IOobject::READ_IF_PRESENT && headerOk())
    )
    {
        // Warn for MUST_READ_IF_MODIFIED
        warnNoRereading<ZoneMesh<ZoneType, MeshType>>();

        PtrList<ZoneType>& zones = *this;

        // Read zones
        Istream& is = readStream(typeName);

        PtrList<entry> patchEntries(is);
        zones.resize(patchEntries.size());

        forAll(zones, zonei)
        {
            zones.set
            (
                zonei,
                ZoneType::New
                (
                    patchEntries[zonei].keyword(),
                    patchEntries[zonei].dict(),
                    zonei,
                    *this
                )
            );
        }

        // Check state of IOstream
        is.check(FUNCTION_NAME);

        close();

        return true;
    }

    // Nothing read
    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ZoneType, class MeshType>
Foam::ZoneMesh<ZoneType, MeshType>::ZoneMesh
(
    const IOobject& io,
    const MeshType& mesh
)
:
    PtrList<ZoneType>(),
    regIOobject(io),
    mesh_(mesh),
    zoneMapPtr_(nullptr)
{
    read();
}


template<class ZoneType, class MeshType>
Foam::ZoneMesh<ZoneType, MeshType>::ZoneMesh
(
    const IOobject& io,
    const MeshType& mesh,
    const label size
)
:
    PtrList<ZoneType>(size),
    regIOobject(io),
    mesh_(mesh),
    zoneMapPtr_(nullptr)
{
    // Optionally read contents, otherwise keep size
    read();
}


template<class ZoneType, class MeshType>
Foam::ZoneMesh<ZoneType, MeshType>::ZoneMesh
(
    const IOobject& io,
    const MeshType& mesh,
    const PtrList<ZoneType>& pzm
)
:
    PtrList<ZoneType>(),
    regIOobject(io),
    mesh_(mesh),
    zoneMapPtr_(nullptr)
{
    if (!read())
    {
        // Nothing read. Use supplied zones
        PtrList<ZoneType>& zones = *this;
        zones.resize(pzm.size());

        forAll(zones, zonei)
        {
            zones.set(zonei, pzm[zonei].clone(*this).ptr());
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class ZoneType, class MeshType>
Foam::ZoneMesh<ZoneType, MeshType>::~ZoneMesh()
{
    clearAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ZoneType, class MeshType>
const Foam::Map<Foam::label>&
Foam::ZoneMesh<ZoneType, MeshType>::zoneMap() const
{
    if (!zoneMapPtr_)
    {
        calcZoneMap();
    }

    return *zoneMapPtr_;
}


template<class ZoneType, class MeshType>
Foam::label Foam::ZoneMesh<ZoneType, MeshType>::whichZone
(
    const label objectIndex
) const
{
    return zoneMap().lookup(objectIndex, -1);
}


template<class ZoneType, class MeshType>
Foam::wordList Foam::ZoneMesh<ZoneType, MeshType>::types() const
{
    return PtrListOps::get<word>(*this, typeOp<ZoneType>());
}


template<class ZoneType, class MeshType>
Foam::wordList Foam::ZoneMesh<ZoneType, MeshType>::names() const
{
    return PtrListOps::get<word>(*this, nameOp<ZoneType>());
}


template<class ZoneType, class MeshType>
Foam::wordList Foam::ZoneMesh<ZoneType, MeshType>::names
(
    const wordRe& matcher
) const
{
    return PtrListOps::names(*this, matcher);
}


template<class ZoneType, class MeshType>
Foam::wordList Foam::ZoneMesh<ZoneType, MeshType>::names
(
    const wordRes& matcher
)
const
{
    return PtrListOps::names(*this, matcher);
}


template<class ZoneType, class MeshType>
Foam::wordList Foam::ZoneMesh<ZoneType, MeshType>::sortedNames() const
{
    wordList sorted(this->names());
    Foam::sort(sorted);

    return sorted;
}


template<class ZoneType, class MeshType>
Foam::wordList Foam::ZoneMesh<ZoneType, MeshType>::sortedNames
(
    const wordRe& matcher
) const
{
    wordList sorted(this->names(matcher));
    Foam::sort(sorted);

    return sorted;
}


template<class ZoneType, class MeshType>
Foam::wordList Foam::ZoneMesh<ZoneType, MeshType>::sortedNames
(
    const wordRes& matcher
)
const
{
    wordList sorted(this->names(matcher));
    Foam::sort(sorted);

    return sorted;
}


template<class ZoneType, class MeshType>
Foam::labelList Foam::ZoneMesh<ZoneType, MeshType>::indices
(
    const wordRe& matcher
) const
{
    if (matcher.empty())
    {
        return labelList();
    }
    return PtrListOps::findMatching(*this, matcher);
}


template<class ZoneType, class MeshType>
Foam::labelList Foam::ZoneMesh<ZoneType, MeshType>::indices
(
    const wordRes& matcher
) const
{
    if (matcher.empty())
    {
        return labelList();
    }
    return PtrListOps::findMatching(*this, matcher);
}


template<class ZoneType, class MeshType>
Foam::label Foam::ZoneMesh<ZoneType, MeshType>::findIndex
(
    const wordRe& key
) const
{
    if (key.empty())
    {
        return -1;
    }
    return PtrListOps::firstMatching(*this, key);
}


template<class ZoneType, class MeshType>
Foam::label Foam::ZoneMesh<ZoneType, MeshType>::findIndex
(
    const wordRes& matcher
) const
{
    if (matcher.empty())
    {
        return -1;
    }
    return PtrListOps::firstMatching(*this, matcher);
}


template<class ZoneType, class MeshType>
Foam::label Foam::ZoneMesh<ZoneType, MeshType>::findZoneID
(
    const word& zoneName
) const
{
    if (zoneName.empty())
    {
        return -1;
    }

    label zoneId = PtrListOps::firstMatching(*this, zoneName);

    if (zoneId < 0)
    {
        DebugInFunction
            << "Zone named " << zoneName << " not found.  "
            << "List of available zone names: " << names() << endl;

        // Used for -dry-run, for example
        if (disallowGenericZones != 0)
        {
            auto& zm = const_cast<ZoneMesh<ZoneType, MeshType>&>(*this);
            zoneId = zm.size();

            Info<< "Creating dummy zone " << zoneName << endl;
            zm.append(new ZoneType(zoneName, zoneId, zm));
        }
    }

    return zoneId;
}


template<class ZoneType, class MeshType>
const ZoneType* Foam::ZoneMesh<ZoneType, MeshType>::cfindZone
(
    const word& zoneName
) const
{
    if (zoneName.empty())
    {
        return nullptr;
    }

    const PtrList<ZoneType>& zones = *this;

    for (auto iter = zones.begin(); iter != zones.end(); ++iter)
    {
        const ZoneType* ptr = iter.get();

        if (ptr && zoneName == ptr->name())
        {
            return ptr;
        }
    }

    // Used for -dry-run, for example
    if (disallowGenericZones != 0)
    {
        auto& zm = const_cast<ZoneMesh<ZoneType, MeshType>&>(*this);

        Info<< "Creating dummy zone " << zoneName << endl;
        zm.append(new ZoneType(zoneName, zm.size(), zm));
    }

    return nullptr;
}


template<class ZoneType, class MeshType>
ZoneType* Foam::ZoneMesh<ZoneType, MeshType>::findZone
(
    const word& zoneName
)
{
    return const_cast<ZoneType*>(this->cfindZone(zoneName));
}


template<class ZoneType, class MeshType>
Foam::bitSet Foam::ZoneMesh<ZoneType, MeshType>::selection
(
    const labelUList& zoneIds
) const
{
    bitSet bitset;

    for (const label zonei : zoneIds)
    {
        #ifdef FULLDEBUG
        if (zonei < 0 || zonei >= this->size())
        {
            FatalErrorInFunction
                << ZoneType::typeName << " "
                << zonei << " out of range [0," << this->size() << ")"
                << abort(FatalError);
        }
        #endif

        bitset.set
        (
            static_cast<const labelList&>(this->operator[](zonei))
        );
    }

    return bitset;
}


template<class ZoneType, class MeshType>
Foam::bitSet Foam::ZoneMesh<ZoneType, MeshType>::selection
(
    const wordRe& matcher
) const
{
    // matcher.empty() is handled by indices()
    return this->selection(this->indices(matcher));
}


template<class ZoneType, class MeshType>
Foam::bitSet Foam::ZoneMesh<ZoneType, MeshType>::selection
(
    const wordRes& matcher
) const
{
    // matcher.empty() is handled by indices()
    return this->selection(this->indices(matcher));
}


template<class ZoneType, class MeshType>
void Foam::ZoneMesh<ZoneType, MeshType>::clearAddressing()
{
    deleteDemandDrivenData(zoneMapPtr_);

    PtrList<ZoneType>& zones = *this;

    for (ZoneType& zn : zones)
    {
        zn.clearAddressing();
    }
}


template<class ZoneType, class MeshType>
void Foam::ZoneMesh<ZoneType, MeshType>::clear()
{
    clearAddressing();
    PtrList<ZoneType>::clear();
}


template<class ZoneType, class MeshType>
bool Foam::ZoneMesh<ZoneType, MeshType>::checkDefinition
(
    const bool report
) const
{
    bool inError = false;

    const PtrList<ZoneType>& zones = *this;

    for (const ZoneType& zn : zones)
    {
        inError |= zn.checkDefinition(report);
    }

    return inError;
}


template<class ZoneType, class MeshType>
bool Foam::ZoneMesh<ZoneType, MeshType>::checkParallelSync
(
    const bool report
) const
{
    if (!Pstream::parRun())
    {
        return false;
    }

    const PtrList<ZoneType>& zones = *this;

    bool hasError = false;

    // Collect all names
    List<wordList> allNames(Pstream::nProcs());
    allNames[Pstream::myProcNo()] = this->names();
    Pstream::gatherList(allNames);
    Pstream::scatterList(allNames);

    List<wordList> allTypes(Pstream::nProcs());
    allTypes[Pstream::myProcNo()] = this->types();
    Pstream::gatherList(allTypes);
    Pstream::scatterList(allTypes);

    // Have every processor check but only master print error.

    for (label proci = 1; proci < allNames.size(); ++proci)
    {
        if
        (
            (allNames[proci] != allNames[0])
         || (allTypes[proci] != allTypes[0])
        )
        {
            hasError = true;

            if (debug || (report && Pstream::master()))
            {
                Info<< " ***Inconsistent zones across processors, "
                       "processor 0 has zone names:" << allNames[0]
                    << " zone types:" << allTypes[0]
                    << " processor " << proci << " has zone names:"
                    << allNames[proci]
                    << " zone types:" << allTypes[proci]
                    << endl;
            }
        }
    }

    // Check contents
    if (!hasError)
    {
        for (const ZoneType& zn : zones)
        {
            if (zn.checkParallelSync(false))
            {
                hasError = true;

                if (debug || (report && Pstream::master()))
                {
                    Info<< " ***Zone " << zn.name()
                        << " of type " << zn.type()
                        << " is not correctly synchronised"
                        << " across coupled boundaries."
                        << " (coupled faces are either not both"
                        << " present in set or have same flipmap)" << endl;
                }
            }
        }
    }

    return hasError;
}


template<class ZoneType, class MeshType>
void Foam::ZoneMesh<ZoneType, MeshType>::movePoints(const pointField& pts)
{
    PtrList<ZoneType>& zones = *this;

    for (ZoneType& zn : zones)
    {
        zn.movePoints(pts);
    }
}


template<class ZoneType, class MeshType>
void Foam::ZoneMesh<ZoneType, MeshType>::updateMetaData()
{
    wordList zoneNames(this->names());
    if (zoneNames.empty())
    {
        this->removeMetaData();
    }
    else
    {
        dictionary& meta = this->getMetaData();
        meta.set("names", zoneNames);
    }
}


template<class ZoneType, class MeshType>
bool Foam::ZoneMesh<ZoneType, MeshType>::writeData(Ostream& os) const
{
    os  << *this;
    return os.good();
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<class ZoneType, class MeshType>
const ZoneType& Foam::ZoneMesh<ZoneType, MeshType>::operator[]
(
    const word& zoneName
) const
{
    const label zonei = findZoneID(zoneName);

    if (zonei < 0)
    {
        FatalErrorInFunction
            << "Zone named " << zoneName << " not found." << nl
            << "Available zone names: " << names() << endl
            << abort(FatalError);
    }

    return operator[](zonei);
}


template<class ZoneType, class MeshType>
ZoneType& Foam::ZoneMesh<ZoneType, MeshType>::operator[]
(
    const word& zoneName
)
{
    const label zonei = findZoneID(zoneName);

    if (zonei < 0)
    {
        FatalErrorInFunction
            << "Zone named " << zoneName << " not found." << nl
            << "Available zone names: " << names() << endl
            << abort(FatalError);
    }

    return operator[](zonei);
}


template<class ZoneType, class MeshType>
ZoneType& Foam::ZoneMesh<ZoneType, MeshType>::operator()
(
    const word& zoneName,
    const bool verbose
)
{
    ZoneType* ptr = findZone(zoneName);

    const bool existing = bool(ptr);

    if (!ptr)
    {
        ptr = new ZoneType(zoneName, this->size(), *this);
        this->append(ptr);
    }

    if (verbose)
    {
        Info<< ZoneType::typeName << ' ' << zoneName
            << " (" << (existing ? "existing" : "new")
            << " at index " << ptr->index() << ')'
            << endl;
    }

    return *ptr;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ZoneType, class MeshType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const ZoneMesh<ZoneType, MeshType>& zones
)
{
    const label sz = zones.size();

    if (sz)
    {
        os  << sz << nl << token::BEGIN_LIST;

        for (label i=0; i < sz; ++i)
        {
            zones[i].writeDict(os);
        }

        os  << token::END_LIST;
    }
    else
    {
        os  << sz << token::BEGIN_LIST << token::END_LIST;
    }

    return os;
}


// ************************************************************************* //
