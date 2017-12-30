/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenCFD Ltd.
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

Description
    read/write auxiliary files for aiding STARCD/OPENFOAM interoperability

    - cellTable information
    - boundaryRegion information
    - interface information

\*----------------------------------------------------------------------------*/

#include "OFstream.H"
#include "ccmReader.H"
#include "ccmInternal.H"  // include last to avoid any strange interactions


// * * * * * * * * * * * * * * * Static Functions  * * * * * * * * * * * * * //

Foam::Map<Foam::word> Foam::ccm::reader::selectPorous
(
    const Map<dictionary>& table
)
{
    Map<word> lookup;

    forAllConstIters(table, iter)
    {
        if (iter().lookupOrDefault<label>("PorosityId", 0) != 0)
        {
            lookup.insert
            (
                iter.key(),
                iter().lookupOrDefault<word>
                (
                    "Label",
                    "cellTable_" + Foam::name(iter.key())
                )
            );
        }
    }

    return lookup;
}


void Foam::ccm::reader::warnDuplicates
(
    const word& context,
    const wordList& lst
)
{
    HashTable<label> hashed(2*lst.size());
    bool duplicates = false;

    for (const word& item : lst)
    {
        // Check duplicate names
        auto iter = hashed.find(item);
        if (iter.found())
        {
            (*iter)++;
            duplicates = true;
        }
        else
        {
            hashed.insert(item, 1);
        }
    }

    // Warn about duplicate names
    if (duplicates)
    {
        Info << nl << "WARNING: " << context << " with identical names:";
        forAllConstIters(hashed, iter)
        {
            if (iter.object() > 1)
            {
                Info << "  " << iter.key();
            }
        }
        Info << nl << endl;
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::ccm::reader::writeInterfaces
(
    const objectRegistry& registry
) const
{
    // Write constant/polyMesh/interface
    IOList<labelList> ioObj
    (
        IOobject
        (
            "interfaces",
            "constant",
            polyMesh::meshSubDir,
            registry,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            false
        )
    );

    ioObj.note() = "as yet unsupported interfaces (baffles)";

    Info<< "Writing " << ioObj.name() << " to " << ioObj.objectPath() << endl;

    OFstream os(ioObj.objectPath());
    ioObj.writeHeader(os);

    os  << bafInterfaces_;
}


// Write List<label> in constant/polyMesh
// - this is crucial for later conversion back to ccm/starcd
void Foam::ccm::reader::writeMeshLabelList
(
    const objectRegistry& registry,
    const word& propertyName,
    const labelList& list,
    IOstream::streamFormat fmt
) const
{
    // Write constant/polyMesh/propertyName
    IOList<label> ioObj
    (
        IOobject
        (
            propertyName,
            "constant",
            polyMesh::meshSubDir,
            registry,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            false
        ),
        list
    );

    ioObj.note() = "persistent data for STARCD <-> OPENFOAM translation";
    Info<< "Writing " << ioObj.name() << " to " << ioObj.objectPath() << endl;

    // NOTE:
    // The cellTableId is an integer and almost always < 1000, thus ASCII
    // will be compacter than binary and makes external scripting easier
    //
    ioObj.writeObject
    (
        fmt,
        IOstream::currentVersion,
        IOstream::UNCOMPRESSED,
        true
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ccm::reader::writeAux
(
    const objectRegistry& registry
) const
{
    cellTable_.writeDict(registry);
    boundaryRegion_.writeDict(registry);
    writeInterfaces(registry);

    // Write origCellId as List<label>
    writeMeshLabelList
    (
        registry,
        "origCellId",
        origCellId_,
        IOstream::BINARY
    );

    // Write cellTableId as List<label>
    // - this is crucial for later conversion back to ccm/starcd
    writeMeshLabelList
    (
        registry,
        "cellTableId",
        cellTableId_,
        IOstream::ASCII
    );
}


// ************************************************************************* //
