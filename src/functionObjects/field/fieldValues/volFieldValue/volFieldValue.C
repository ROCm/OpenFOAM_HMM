/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

#include "volFieldValue.H"
#include "fvMesh.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace fieldValues
{
    defineTypeNameAndDebug(volFieldValue, 0);
    addToRunTimeSelectionTable(fieldValue, volFieldValue, dictionary);
    addToRunTimeSelectionTable(functionObject, volFieldValue, dictionary);
}
}
}

const Foam::Enum
<
    Foam::functionObjects::fieldValues::volFieldValue::operationType
>
Foam::functionObjects::fieldValues::volFieldValue::operationTypeNames_
{
    // Normal operations
    { operationType::opNone, "none" },
    { operationType::opMin, "min" },
    { operationType::opMax, "max" },
    { operationType::opSum, "sum" },
    { operationType::opSumMag, "sumMag" },
    { operationType::opAverage, "average" },
    { operationType::opVolAverage, "volAverage" },
    { operationType::opVolIntegrate, "volIntegrate" },
    { operationType::opCoV, "CoV" },

    // Using weighting
    { operationType::opWeightedSum, "weightedSum" },
    { operationType::opWeightedAverage, "weightedAverage" },
    { operationType::opWeightedVolAverage, "weightedVolAverage" },
    { operationType::opWeightedVolIntegrate, "weightedVolIntegrate" },
};


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::fieldValues::volFieldValue::usesVol() const
{
    // Only a few operations require the cell volume
    switch (operation_)
    {
        case opVolAverage:
        case opVolIntegrate:
        case opWeightedVolAverage:
        case opWeightedVolIntegrate:
        case opCoV:
            return true;

        default:
            return false;
    }
}


bool Foam::functionObjects::fieldValues::volFieldValue::usesWeight() const
{
    // Operation specifically tagged to require a weight field
    return (operation_ & typeWeighted);
}


bool Foam::functionObjects::fieldValues::volFieldValue::canWeight
(
    const scalarField& weightField
) const
{
    return
    (
        usesWeight()
     && returnReduce(!weightField.empty(), orOp<bool>()) // On some processor
    );
}


void Foam::functionObjects::fieldValues::volFieldValue::initialise
(
    const dictionary& dict
)
{
    weightFieldName_ = "none";
    if (usesWeight())
    {
        if (dict.readIfPresent("weightField", weightFieldName_))
        {
            Info<< "    weight field = " << weightFieldName_;
        }
        else
        {
            // Suggest possible alternative unweighted operation?
            FatalIOErrorInFunction(dict)
                << "The '" << operationTypeNames_[operation_]
                << "' operation is missing a weightField." << nl
                << "Either provide the weightField, "
                << "use weightField 'none' to suppress weighting," << nl
                << "or use a different operation."
                << exit(FatalIOError);
        }
    }

    Info<< nl << endl;
}


void Foam::functionObjects::fieldValues::volFieldValue::writeFileHeader
(
    Ostream& os
) const
{
    volRegion::writeFileHeader(*this, os);
    if (weightFieldName_ != "none")
    {
        writeHeaderValue(os, "Weight field", weightFieldName_);
    }

    writeCommented(os, "Time");

    for (const word& fieldName : fields_)
    {
        os  << tab << operationTypeNames_[operation_]
            << "(" << fieldName << ")";
    }

    os  << endl;
}


Foam::label Foam::functionObjects::fieldValues::volFieldValue::writeAll
(
    const scalarField& V,
    const scalarField& weightField
)
{
    label nProcessed = 0;

    for (const word& fieldName : fields_)
    {
        if
        (
            writeValues<scalar>(fieldName, V, weightField)
         || writeValues<vector>(fieldName, V, weightField)
         || writeValues<sphericalTensor>(fieldName, V, weightField)
         || writeValues<symmTensor>(fieldName, V, weightField)
         || writeValues<tensor>(fieldName, V, weightField)
        )
        {
            ++nProcessed;
        }
        else
        {
            WarningInFunction
                << "Requested field " << fieldName
                << " not found in database and not processed"
                << endl;
        }
    }

    return nProcessed;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldValues::volFieldValue::volFieldValue
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldValue(name, runTime, dict, typeName),
    volRegion(fieldValue::mesh_, dict),
    operation_(operationTypeNames_.lookup("operation", dict)),
    weightFieldName_("none")
{
    read(dict);
    writeFileHeader(file());
}


Foam::functionObjects::fieldValues::volFieldValue::volFieldValue
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    fieldValue(name, obr, dict, typeName),
    volRegion(fieldValue::mesh_, dict),
    operation_(operationTypeNames_.lookup("operation", dict)),
    weightFieldName_("none")
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::fieldValues::volFieldValue::~volFieldValue()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldValues::volFieldValue::read
(
    const dictionary& dict
)
{
    fieldValue::read(dict);
    initialise(dict);

    return true;
}


bool Foam::functionObjects::fieldValues::volFieldValue::write()
{
    fieldValue::write();

    if (Pstream::master())
    {
        writeTime(file());
    }

    // Only some operations need the cell volume
    scalarField V;
    if (usesVol())
    {
        V = filterField(fieldValue::mesh_.V());
    }

    // Weight field - zero-size means weight = 1
    scalarField weightField;
    if (weightFieldName_ != "none")
    {
        weightField = getFieldValues<scalar>(weightFieldName_, true);
    }

    writeAll(V, weightField);

    if (Pstream::master())
    {
        file()<< endl;
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
