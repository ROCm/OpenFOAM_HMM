/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd
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

#include "fieldValueDelta.H"
#include "ListOps.H"
#include "Time.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace fieldValues
    {
        defineTypeNameAndDebug(fieldValueDelta, 0);
    }

    template<>
    const char*
    NamedEnum<fieldValues::fieldValueDelta::operationType, 5>::names[] =
    {
        "add",
        "subtract",
        "min",
        "max",
        "average"
    };

    const NamedEnum<fieldValues::fieldValueDelta::operationType, 5>
        fieldValues::fieldValueDelta::operationTypeNames_;
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::fieldValues::fieldValueDelta::writeFileHeader(Ostream& os) const
{
    const wordList& fields1 = source1Ptr_->fields();
    const wordList& fields2 = source2Ptr_->fields();

    DynamicList<word> commonFields(fields1.size());
    forAll(fields1, i)
    {
        label index = findIndex(fields2, fields1[i]);
        if (index != -1)
        {
            commonFields.append(fields1[i]);
        }
    }

    writeHeaderValue(os, "Source1", source1Ptr_->name());
    writeHeaderValue(os, "Source2", source2Ptr_->name());
    writeHeaderValue(os, "Operation", operationTypeNames_[operation_]);
    writeCommented(os, "Time");

    forAll(commonFields, i)
    {
        os  << tab << commonFields[i];
    }

    os  << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fieldValues::fieldValueDelta::fieldValueDelta
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    functionObjectState(obr, name),
    functionObjectFile(obr, name, typeName, dict),
    obr_(obr),
    loadFromFiles_(loadFromFiles),
    log_(true),
    operation_(opSubtract),
    source1Ptr_(NULL),
    source2Ptr_(NULL)
{
    if (setActive<fvMesh>())
    {
        read(dict);
        writeFileHeader(file());
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fieldValues::fieldValueDelta::~fieldValueDelta()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fieldValues::fieldValueDelta::read(const dictionary& dict)
{
    if (active_)
    {
        functionObjectFile::read(dict);

        log_ = dict.lookupOrDefault<Switch>("log", true);
        source1Ptr_.reset
        (
            fieldValue::New
            (
                name_ + ".source1",
                obr_,
                dict.subDict("source1"),
                loadFromFiles_,
                false
            ).ptr()
        );
        source2Ptr_.reset
        (
            fieldValue::New
            (
                name_ + ".source2",
                obr_,
                dict.subDict("source2"),
                loadFromFiles_,
                false
            ).ptr()
        );

        operation_ = operationTypeNames_.read(dict.lookup("operation"));
    }
}


void Foam::fieldValues::fieldValueDelta::write()
{
    // Do nothing
}


void Foam::fieldValues::fieldValueDelta::execute()
{
    if (active_)
    {
        source1Ptr_->write();
        source2Ptr_->write();

        writeTime(file());

        if (log_) Info << type() << " " << name_ << " output:" << endl;

        const word& name1 = source1Ptr_->name();
        const word& name2 = source2Ptr_->name();

        const wordList entries1 = objectResultEntries(name1);
        const wordList entries2 = objectResultEntries(name2);

        if (entries1.size() != entries2.size())
        {
            FatalErrorInFunction
                << name_ << ": objects must generate the same number of results"
                << nl
                << "    " << name1 << " objects: " << entries1 << nl
                << "    " << name2 << " objects: " << entries2 << nl
                << exit(FatalError);
        }

        forAll(entries1, i)
        {
            const word& entry1(entries1[i]);
            const word& entry2(entries2[i]);
            const word type1 = objectResultType(name1, entry1);
            const word type2 = objectResultType(name2, entry2);

            if (type1 != type2)
            {
                FatalErrorInFunction
                    << name_
                    << ": input values for operation must be of the same type"
                    << nl
                    << "    " << entry1 << ": " << type1 << nl
                    << "    " << entry2 << ": " << type2 << nl
                    << exit(FatalError);
            }

            bool found = false;

            apply<scalar>(type1, name1, name2, entry1, entry2, found);
            apply<vector>(type1, name1, name2, entry1, entry2, found);
            apply<sphericalTensor>(type1, name1, name2, entry1, entry2, found);
            apply<symmTensor>(type1, name1, name2, entry1, entry2, found);
            apply<tensor>(type1, name1, name2, entry1, entry2, found);

            if (log_ && !found)
            {
                Info<< "Operation between "
                    << name1 << " with result " << entry1 << " and "
                    << name2 << " with result " << entry2 << " not applied"
                    << endl;
            }
        }

        if (log_)
        {
            if (entries1.empty())
            {
                Info<< "    none";
            }
            Info<< endl;
        }

        file()<< endl;
    }
}


void Foam::fieldValues::fieldValueDelta::end()
{
    if (active_)
    {
        execute();
    }
}


void Foam::fieldValues::fieldValueDelta::timeSet()
{
    // Do nothing
}


void Foam::fieldValues::fieldValueDelta::updateMesh(const mapPolyMesh&)
{
    // Do nothing
}


void Foam::fieldValues::fieldValueDelta::movePoints(const polyMesh&)
{
    // Do nothing
}


// ************************************************************************* //
