/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "cellTable.H"
#include "IOMap.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::cellTable::defaultMaterial_ = "fluid";


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::Map<Foam::label> Foam::cellTable::zoneMap() const
{
    Map<label> lookup;

    label zonei = 0;
    forAllConstIters(*this, iter)
    {
        lookup.insert(iter.key(), zonei++);
    }

    return lookup;
}


Foam::wordList Foam::cellTable::namesList() const
{
    Map<word> lookup = names();
    wordList list(lookup.size());

    label zonei = 0;
    forAllConstIters(lookup, iter)
    {
        list[zonei++] = *iter;
    }

    return list;
}


void Foam::cellTable::addDefaults()
{
    forAllIters(*this, iter)
    {
        if (!iter().found("MaterialType"))
        {
            iter().add("MaterialType", defaultMaterial_);
        }
    }
}


void Foam::cellTable::setEntry
(
    const label id,
    const word& keyWord,
    const word& value
)
{
    dictionary dict;
    dict.add(keyWord, value);

    iterator iter = find(id);
    if (iter.found())
    {
        iter().merge(dict);
    }
    else
    {
        insert(id, dict);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellTable::cellTable()
:
    Map<dictionary>()
{}


Foam::cellTable::cellTable
(
    const objectRegistry& registry,
    const word& name,
    const fileName& instance
)
:
    Map<dictionary>()
{
    readDict(registry, name, instance);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::label Foam::cellTable::append(const dictionary& dict)
{
    label maxId = -1;
    forAllConstIters(*this, iter)
    {
        if (maxId < iter.key())
        {
            maxId = iter.key();
        }
    }

    insert(++maxId, dict);
    return maxId;
}


Foam::Map<Foam::word> Foam::cellTable::names() const
{
    Map<word> lookup;

    forAllConstIters(*this, iter)
    {
        lookup.insert
        (
            iter.key(),
            iter().getOrDefault<word>
            (
                "Label",
                "cellTable_" + Foam::name(iter.key())
            )
        );
    }

    return lookup;
}


Foam::Map<Foam::word> Foam::cellTable::names
(
    const wordRes& patterns
) const
{
    Map<word> lookup;

    forAllConstIters(*this, iter)
    {
        const word lookupName = iter().getOrDefault<word>
        (
            "Label",
            "cellTable_" + Foam::name(iter.key())
        );

        if (patterns.match(lookupName))
        {
            lookup.insert(iter.key(), lookupName);
        }
    }

    return lookup;
}


Foam::word Foam::cellTable::name(const label id) const
{
    word theName("cellTable_" + Foam::name(id));

    const_iterator iter = cfind(id);
    if (iter.found())
    {
        iter().readIfPresent("Label", theName);
    }

    return theName;
}


Foam::label Foam::cellTable::findIndex(const word& name) const
{
    if (name.empty())
    {
        return -1;
    }

    forAllConstIters(*this, iter)
    {
        if (iter().getOrDefault<word>("Label", word::null) == name)
        {
            return iter.key();
        }
    }

    return -1;
}


Foam::Map<Foam::word> Foam::cellTable::materialTypes() const
{
    Map<word> lookup;

    forAllConstIters(*this, iter)
    {
        lookup.insert
        (
            iter.key(),
            iter().getOrDefault<word>("MaterialType", defaultMaterial_)
        );
    }

    return lookup;
}


Foam::Map<Foam::word> Foam::cellTable::selectType(const word& matl) const
{
    Map<word> lookup;

    forAllConstIters(*this, iter)
    {
        const label index = iter.key();
        const dictionary& dict = iter.val();

        if
        (
            matl
         == dict.getOrDefault<word>("MaterialType", defaultMaterial_)
        )
        {
            lookup.insert
            (
                index,
                dict.getOrDefault<word>
                (
                    "Label",
                    "cellTable_" + Foam::name(iter.key())
                )
            );
        }
    }

    return lookup;
}


Foam::Map<Foam::word> Foam::cellTable::fluids() const
{
    return selectType("fluid");
}


Foam::Map<Foam::word> Foam::cellTable::solids() const
{
    return selectType("solid");
}


Foam::Map<Foam::word> Foam::cellTable::shells() const
{
    return selectType("shell");
}



void Foam::cellTable::setMaterial(const label id, const word& matlType)
{
    setEntry(id, "MaterialType", matlType);
}


void Foam::cellTable::setName(const label id, const word& name)
{
    setEntry(id, "Label", name);
}


void Foam::cellTable::setName(const label id)
{
    iterator iter = find(id);

    if (!iter.found() || !iter().found("Label"))
    {
        setName(id, "cellTable_" + Foam::name(id));
    }
}


void Foam::cellTable::readDict
(
    const objectRegistry& registry,
    const word& name,
    const fileName& instance
)
{
    clear();

    // read constant/dictName
    IOMap<dictionary> ioObj
    (
        IOobject
        (
            name,
            instance,
            registry,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        )
    );

    if (ioObj.headerOk())
    {
        *this = ioObj;
        addDefaults();
    }
    else
    {
        Info<< "no constant/cellTable information available" << endl;
    }
}


void Foam::cellTable::writeDict
(
    const objectRegistry& registry,
    const word& name,
    const fileName& instance
) const
{
    // write constant/dictName
    IOMap<dictionary> ioObj
    (
        IOobject
        (
            name,
            instance,
            registry,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    ioObj.note() =
        "persistent data for thirdParty mesh <-> OpenFOAM translation";

    Info<< "Writing " << ioObj.name() << " to " << ioObj.objectPath() << endl;

    OFstream os(ioObj.objectPath());
    ioObj.writeHeader(os);
    os << *this;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::cellTable::operator=(const cellTable& rhs)
{
    Map<dictionary>::operator=(rhs);
    addDefaults();
}


void Foam::cellTable::operator=(const Map<dictionary>& rhs)
{
    Map<dictionary>::operator=(rhs);
    addDefaults();
}


void Foam::cellTable::operator=(const polyMesh& mesh)
{
    Map<dictionary> zoneDict;

    // create cellTableId and cellTable based on cellZones
    label nZoneCells = 0;

    wordList zoneNames = mesh.cellZones().names();
    label unZonedType = zoneNames.size() + 1;

    // do cell zones
    forAll(mesh.cellZones(), zoneI)
    {
        const cellZone& cZone = mesh.cellZones()[zoneI];
        nZoneCells += cZone.size();

        dictionary dict;
        dict.add("Label", zoneNames[zoneI]);
        zoneDict.insert(zoneI + 1, dict);
    }

    // collect unzoned cells
    // special case: no zones at all - do entire mesh
    if (nZoneCells == 0)
    {
        zoneDict.clear();
        unZonedType = 1;
    }

    if (mesh.nCells() > nZoneCells)
    {
        zoneDict.insert
        (
            unZonedType,
            dictionary(IStringStream("Label cells;")())
        );
    }

    Map<dictionary>::operator=(zoneDict);
    addDefaults();
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

void Foam::cellTable::addCellZones
(
    polyMesh& mesh,
    const labelList& tableIds
) const
{
    Map<label> typeToZone = zoneMap();
    List<DynamicList<label>> zoneCells(size());

    forAll(tableIds, celli)
    {
        const auto iter = typeToZone.cfind(tableIds[celli]);
        if (iter.found())
        {
            zoneCells[*iter].append(celli);
        }
    }

    // track which zones were actually used
    labelList zoneUsed(zoneCells.size());
    wordList  zoneNames(namesList());

    label nZone = 0;
    forAll(zoneCells, zoneI)
    {
        zoneCells[zoneI].shrink();
        if (zoneCells[zoneI].size())
        {
            zoneUsed[nZone++] = zoneI;
        }
    }
    zoneUsed.setSize(nZone);

    cellZoneMesh& czMesh = mesh.cellZones();

    czMesh.clear();
    if (nZone <= 1)
    {
        Info<< "cellZones not used" << endl;
        return;
    }
    czMesh.setSize(nZone);

    forAll(zoneUsed, zoneI)
    {
        const label origZoneI = zoneUsed[zoneI];

        Info<< "cellZone " << zoneI
            << " (size: "  << zoneCells[origZoneI].size()
            << ") name: "  << zoneNames[origZoneI] << endl;

        czMesh.set
        (
            zoneI,
            new cellZone
            (
                zoneNames[origZoneI],
                zoneCells[origZoneI],
                zoneI,
                czMesh
            )
        );
    }
    czMesh.writeOpt(IOobject::AUTO_WRITE);
}


void Foam::cellTable::combine(const dictionary& mapDict, labelList& tableIds)
{
    if (mapDict.empty())
    {
        return;
    }

    Map<word> origNames(names());
    labelList mapping(identity(max(origNames.toc()) + 1));

    bool remap = false;
    forAllConstIters(mapDict, iter)
    {
        wordRes patterns(iter().stream());

        // find all matches
        Map<word> matches;
        forAllConstIters(origNames, namesIter)
        {
            if (patterns.match(namesIter()))
            {
                matches.insert(namesIter.key(), namesIter());
            }
        }

        if (matches.size())
        {
            label targetId = this->findIndex(iter().keyword());

            Info<< "combine cellTable: " << iter().keyword();
            if (targetId < 0)
            {
                // not found - reuse 1st element but with different name
                targetId = min(matches.toc());
                operator[](targetId).set("Label", iter().keyword());

                Info<< " = (";
            }
            else
            {
                Info<< " += (";
            }


            // the mapping and name for targetId is already okay
            matches.erase(targetId);
            origNames.erase(targetId);

            // remove matched names, leaving targetId on 'this'
            this->erase(matches);
            origNames.erase(matches);

            forAllConstIters(matches, matchIter)
            {
                mapping[matchIter.key()] = targetId;
                Info<< " " << matchIter();
            }
            Info<< " )" << endl;

            remap = true;
        }
    }

    if (remap)
    {
        inplaceRenumber(mapping, tableIds);
    }
}

// ************************************************************************* //
