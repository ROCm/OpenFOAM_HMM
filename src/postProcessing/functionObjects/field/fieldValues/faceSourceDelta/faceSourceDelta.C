/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "faceSourceDelta.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::fieldValues::faceSourceDelta, 0);

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fieldValues::faceSourceDelta::updateMesh(const mapPolyMesh&)
{
    // Do nothing
}


void Foam::fieldValues::faceSourceDelta::movePoints(const Field<point>&)
{
    // Do nothing
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fieldValues::faceSourceDelta::faceSourceDelta
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    functionObjectFile(obr, name, typeName),
    obr_(obr),
    log_(false),
    faceSource1_
    (
        name + ".faceSource1",
        obr,
        dict.subDict("faceSource1"),
        loadFromFiles
    ),
    faceSource2_
    (
        name + ".faceSource2",
        obr,
        dict.subDict("faceSource2"),
        loadFromFiles
    )
{
    read(dict);
}


void Foam::fieldValues::faceSourceDelta::writeFileHeader(const label i)
{
    const wordList& fields1 = faceSource1_.fields();
    const wordList& fields2 = faceSource2_.fields();

    DynamicList<word> commonFields(fields1.size());
    forAll(fields1, i)
    {
        label index = findIndex(fields2, fields1[i]);
        if (index != -1)
        {
            commonFields.append(fields1[i]);
        }
    }

    file() << "# Time";

    forAll(commonFields, i)
    {
        file()<< tab << commonFields[i];
    }

    file() << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fieldValues::faceSourceDelta::~faceSourceDelta()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fieldValues::faceSourceDelta::read(const dictionary& dict)
{
    log_ = dict.lookupOrDefault<Switch>("log", false);
    faceSource1_.read(dict.subDict("faceSource1"));
    faceSource2_.read(dict.subDict("faceSource2"));
}


void Foam::fieldValues::faceSourceDelta::write()
{
    functionObjectFile::write();

    faceSource1_.write();
    faceSource2_.write();

    if (Pstream::master())
    {
        file()<< obr_.time().timeName();
    }

    if (log_)
    {
        Info<< type() << " output:" << endl;
    }

    bool found = false;
    processFields<scalar>(found);
    processFields<vector>(found);
    processFields<sphericalTensor>(found);
    processFields<symmTensor>(found);
    processFields<tensor>(found);

    if (Pstream::master())
    {
        file()<< endl;
    }

    if (log_)
    {
        if (!found)
        {
            Info<< "    none" << endl;
        }
        else
        {
            Info<< endl;
        }
    }
}


void Foam::fieldValues::faceSourceDelta::execute()
{
    // Do nothing
}


void Foam::fieldValues::faceSourceDelta::end()
{
    // Do nothing
}


// ************************************************************************* //
