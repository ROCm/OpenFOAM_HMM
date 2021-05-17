/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
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

#include "multiFieldValue.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace fieldValues
{
    defineTypeNameAndDebug(multiFieldValue, 0);
    addToRunTimeSelectionTable(functionObject, multiFieldValue, dictionary);
}
}
}


const Foam::Enum
<
    Foam::functionObjects::fieldValues::multiFieldValue::operationType
>
Foam::functionObjects::fieldValues::multiFieldValue::operationTypeNames_
({
    { operationType::opSum, "sum" },
    { operationType::opAdd, "add" },
    { operationType::opSubtract, "subtract" },
    { operationType::opMin, "min" },
    { operationType::opMax, "max" },
    { operationType::opAverage, "average" },
});


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::fieldValues::multiFieldValue::writeFileHeader
(
    Ostream& os
) const
{
    const wordList& fields0 = functions_[0].fields();

    DynamicList<word> commonFields(fields0.size());

    for (const word& fieldName : fields0)
    {
        bool common = true;

        for (label functioni=1; functioni < functions_.size(); ++functioni)
        {
            if (!functions_[functioni].fields().found(fieldName))
            {
                common = false;
                break;
            }
        }

        if (common)
        {
            commonFields.append(fieldName);
        }
    }

    forAll(functions_, functioni)
    {
        writeHeaderValue
        (
            os,
            "Source" + Foam::name(functioni),
            functions_[functioni].name()
        );
    }

    writeHeaderValue(os, "Operation", operationTypeNames_[operation_]);
    writeCommented(os, "Time");

    for (const word& fieldName : commonFields)
    {
        os  << tab << fieldName;
    }

    os  << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldValues::multiFieldValue::multiFieldValue
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    stateFunctionObject(name, runTime),
    writeFile(runTime, name, typeName, dict),
    operation_(opSubtract),
    functions_()
{
    if (read(dict))
    {
        writeFileHeader(file());
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldValues::multiFieldValue::read
(
    const dictionary& dict
)
{
    if (stateFunctionObject::read(dict) && writeFile::read(dict))
    {
        const dictionary& functionsDict = dict.subDict("functions");
        functions_.resize(functionsDict.size());

        if (functions_.empty())
        {
            WarningInFunction
                << "No functions specified"
                << endl;
            return false;
        }

        label functioni = 0;
        for (const entry& dEntry : functionsDict)
        {
            if (!dEntry.isDict())
            {
                FatalIOErrorInFunction(dict)
                    << "Functions must be specified in dictionary format"
                    << exit(FatalIOError);
            }

            const dictionary& localDict = dEntry.dict();

            functions_.set
            (
                functioni,
                fieldValue::New
                (
                    IOobject::scopedName(name(), localDict.dictName()),
                    time(),
                    localDict,
                    false
                )
            );

            ++functioni;
        }

        operation_ = operationTypeNames_.get("operation", dict);

        return true;
    }

    return false;
}


bool Foam::functionObjects::fieldValues::multiFieldValue::write()
{
    if (functions_.empty())
    {
        return false;
    }

    Log << type() << " " << name() << " write:" << endl;

    const label nFunction = functions_.size();
    wordList entries0;
    label nEntries = -1;

    wordList names(nFunction);
    List<wordList> entries;
    List<wordList> types;

    forAll(functions_, functioni)
    {
        auto& f = functions_[functioni];
        names[functioni] = f.name();

        // Note: results are not available until the call to write()
        f.write();

        const wordList e(objectResultEntries(f.name()));

        if (functioni == 0)
        {
            entries0 = e;
            nEntries = e.size();
            entries.resize(nEntries);
            types.resize(nEntries);

            forAll(entries, entryi)
            {
                entries[entryi].resize(nFunction);
                types[entryi].resize(nFunction);
            }
        }

        if (e.size() != nEntries)
        {
            const word& f0Name = functions_[0].name();

            FatalErrorInFunction
                << "Inconsistent number of result entries" << nl
                << "    " << f0Name << " entries:" << entries0 << nl
                << "    " << f.name() << " entries:" << e
                << abort(FatalError);
        }

        forAll(e, entryi)
        {
            entries[entryi][functioni] = e[entryi];
            types[entryi][functioni] = objectResultType(f.name(), e[entryi]);
        }
    }

    writeCurrentTime(file());

    forAll(entries, entryi)
    {
        const wordList& entriesi = entries[entryi];
        const word& t0 = types[entryi][0];
        const wordList& typesi = types[entryi];
        forAll(typesi, functioni)
        {
            const word& t = typesi[functioni];

            if (t != t0)
            {
                FatalErrorInFunction
                    << "Inconsistent function result types" << nl
                    << "    " << functions_[0].name()
                    << " result type:" << t0 << nl
                    << "    " << functions_[functioni].name()
                    << " result type:" << typesi[functioni]
                    << abort(FatalError);
            }
        }

        const bool ok
        (
            applyOperation<scalar>(t0, names, entriesi)
         || applyOperation<vector>(t0, names, entriesi)
         || applyOperation<sphericalTensor>(t0, names, entriesi)
         || applyOperation<symmTensor>(t0, names, entriesi)
         || applyOperation<tensor>(t0, names, entriesi)
        );

        if (!ok)
        {
            Log << "Operation not applied between functions:" << nl
                << flatOutput(names, FlatOutput::BareComma{}) << nl
                << "with result names:" << nl
                << flatOutput(entriesi, FlatOutput::BareComma{})
                << endl;
        }
    }

    Log << (nEntries == 0 ? "    none" : "") << endl;

    file()<< endl;

    return true;
}


bool Foam::functionObjects::fieldValues::multiFieldValue::execute()
{
    return true;
}


// ************************************************************************* //
