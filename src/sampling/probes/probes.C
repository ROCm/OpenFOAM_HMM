/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "probes.H"
#include "volFields.H"
#include "dictionary.H"
#include "Time.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::probes, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::probes::findCells(const fvMesh& mesh)
{
    cellList_.clear();
    cellList_.setSize(size());

    forAll(*this, probeI)
    {
        const vector& location = operator[](probeI);

        cellList_[probeI] = mesh.findCell(location);

        if (debug && cellList_[probeI] != -1)
        {
            Pout<< "probes : found point " << location
                << " in cell " << cellList_[probeI] << endl;
        }
    }


    // Check if all probes have been found.
    forAll(cellList_, probeI)
    {
        const vector& location = operator[](probeI);
        label cellI = cellList_[probeI];

        // Check at least one processor with cell.
        reduce(cellI, maxOp<label>());

        if (cellI == -1)
        {
            if (Pstream::master())
            {
                WarningIn("probes::read()")
                    << "Did not find location " << location
                    << " in any cell. Skipping location." << endl;
            }
        }
        else
        {
            // Make sure location not on two domains.
            if (cellList_[probeI] != -1 && cellList_[probeI] != cellI)
            {
                WarningIn("probes::read()")
                    << "Location " << location
                    << " seems to be on multiple domains:"
                    << " cell " << cellList_[probeI]
                    << " on my domain " << Pstream::myProcNo()
                        << " and cell " << cellI << " on some other domain."
                    << endl
                    << "This might happen if the probe location is on"
                    << " a processor patch. Change the location slightly"
                    << " to prevent this." << endl;
            }
        }
    }
}


Foam::label Foam::probes::prepare()
{
    const label nFields = classifyFields();

    // adjust file streams
    if (Pstream::master())
    {
        wordHashSet currentFields;

        forAll(scalarFields_, fieldI)
        {
            currentFields.set(scalarFields_[fieldI]);
        }
        forAll(vectorFields_, fieldI)
        {
            currentFields.set(vectorFields_[fieldI]);
        }
        forAll(sphericalTensorFields_, fieldI)
        {
            currentFields.set(sphericalTensorFields_[fieldI]);
        }
        forAll(symmTensorFields_, fieldI)
        {
            currentFields.set(symmTensorFields_[fieldI]);
        }
        forAll(tensorFields_, fieldI)
        {
            currentFields.set(tensorFields_[fieldI]);
        }

        if (debug)
        {
            Info<< "Probing fields:" << currentFields << nl
                << "Probing locations:" << *this << nl
                << endl;
        }


        fileName probeDir;
        fileName probeSubDir = name_;

        if (mesh_.name() != polyMesh::defaultRegion)
        {
            probeSubDir = probeSubDir/mesh_.name();
        }
        probeSubDir = probeSubDir/mesh_.time().timeName();

        if (Pstream::parRun())
        {
            // Put in undecomposed case
            // (Note: gives problems for distributed data running)
            probeDir = mesh_.time().path()/".."/probeSubDir;
        }
        else
        {
            probeDir = mesh_.time().path()/probeSubDir;
        }

        // ignore known fields, close streams for fields that no longer exist
        forAllIter(HashPtrTable<OFstream>, probeFilePtrs_, iter)
        {
            if (!currentFields.erase(iter.key()))
            {
                if (debug)
                {
                    Info<< "close probe stream: " << iter()->name() << endl;
                }

                delete probeFilePtrs_.remove(iter);
            }
        }

        // currentFields now just has the new fields - open streams for them
        forAllConstIter(wordHashSet, currentFields, iter)
        {
            const word& fieldName = iter.key();

            // Create directory if does not exist.
            mkDir(probeDir);

            OFstream* sPtr = new OFstream(probeDir/fieldName);

            if (debug)
            {
                Info<< "open  probe stream: " << sPtr->name() << endl;
            }

            probeFilePtrs_.insert(fieldName, sPtr);

            unsigned int w = IOstream::defaultPrecision() + 7;

            for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
            {
                *sPtr<< '#' << setw(IOstream::defaultPrecision() + 6)
                    << vector::componentNames[cmpt];

                forAll(*this, probeI)
                {
                    *sPtr<< ' ' << setw(w) << operator[](probeI)[cmpt];
                }
                *sPtr << endl;
            }

            *sPtr<< '#' << setw(IOstream::defaultPrecision() + 6)
                << "Time" << endl;
        }
    }

    return nFields;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::probes::probes
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    pointField(0),
    name_(name),
    mesh_(refCast<const fvMesh>(obr)),
    loadFromFiles_(loadFromFiles)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::probes::~probes()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::probes::execute()
{
    // Do nothing - only valid on write
}


void Foam::probes::end()
{
    // Do nothing - only valid on write
}


void Foam::probes::write()
{
    if (size() && prepare())
    {
        sampleAndWrite(scalarFields_);
        sampleAndWrite(vectorFields_);
        sampleAndWrite(sphericalTensorFields_);
        sampleAndWrite(symmTensorFields_);
        sampleAndWrite(tensorFields_);
    }
}


void Foam::probes::read(const dictionary& dict)
{
    dict.lookup("probeLocations") >> *this;
    dict.lookup("fields") >> fieldSelection_;

    // redetermined all cell locations
    findCells(mesh_);
    prepare();
}


// ************************************************************************* //
