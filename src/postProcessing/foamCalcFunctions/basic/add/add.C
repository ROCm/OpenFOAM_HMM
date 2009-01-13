/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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

#include "add.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace calcTypes
    {
        defineTypeNameAndDebug(add, 0);
        addToRunTimeSelectionTable(calcType, add, dictionary);
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::calcTypes::add::writeAddFields
(
    const Time& runTime,
    const fvMesh& mesh,
    const IOobject& baseFieldHeader
)
{
    bool processed = false;

    IOobject addFieldHeader
    (
        addFieldName_,
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    if (addFieldHeader.headerOk())
    {
        writeAddField<scalar>
        (
            baseFieldHeader,
            addFieldHeader,
            mesh,
            processed
        );
        writeAddField<vector>
        (
            baseFieldHeader,
            addFieldHeader,
            mesh,
            processed
        );
        writeAddField<sphericalTensor>
        (
            baseFieldHeader,
            addFieldHeader,
            mesh,
            processed
        );
        writeAddField<symmTensor>
        (
            baseFieldHeader,
            addFieldHeader,
            mesh,
            processed
        );
        writeAddField<tensor>
        (
            baseFieldHeader,
            addFieldHeader,
            mesh,
            processed
        );

        if (!processed)
        {
            FatalError
                << "Unable to process " << baseFieldName_
                << " + " << addFieldName_ << nl
                << "No call to add for fields of type "
                << baseFieldHeader.headerClassName() << " + "
                << addFieldHeader.headerClassName() << nl << nl
                << exit(FatalError);
        }
    }
    else
    {
        FatalErrorIn("calcTypes::add::writeAddFields()")
            << "Unable to read add field: " << addFieldName_
            << nl << exit(FatalError);
    }
}


void Foam::calcTypes::add::writeAddValues
(
    const Time& runTime,
    const fvMesh& mesh,
    const IOobject& baseFieldHeader
)
{
    bool processed = false;

    writeAddValue<scalar>
    (
        baseFieldHeader,
        addValueStr_,
        mesh,
        processed
    );
    writeAddValue<vector>
    (
        baseFieldHeader,
        addValueStr_,
        mesh,
        processed
    );
    writeAddValue<sphericalTensor>
    (
        baseFieldHeader,
        addValueStr_,
        mesh,
        processed
    );
    writeAddValue<symmTensor>
    (
        baseFieldHeader,
        addValueStr_,
        mesh,
        processed
    );
    writeAddValue<tensor>
    (
        baseFieldHeader,
        addValueStr_,
        mesh,
        processed
    );

    if (!processed)
    {
        FatalErrorIn("calcTypes::add::writeAddValues()")
            << "Unable to process " << baseFieldName_
            << " + " << addValueStr_ << nl
            << "No call to add for fields of type "
            << baseFieldHeader.headerClassName() << nl << nl
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::calcTypes::add::add()
:
    calcType(),
    baseFieldName_(""),
    calcType_(FIELD),
    addFieldName_(""),
    addValueStr_(""),
    resultName_("")
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::calcTypes::add::~add()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::calcTypes::add::init()
{
    argList::validArgs.append("add");
    argList::validArgs.append("baseField");
    argList::validOptions.insert("field", "fieldName");
    argList::validOptions.insert("value", "valueString");
    argList::validOptions.insert("resultName", "fieldName");
}


void Foam::calcTypes::add::preCalc
(
    const argList& args,
    const Time& runTime,
    const fvMesh& mesh
)
{
    baseFieldName_ = args.additionalArgs()[1];

    if (args.options().found("field"))
    {
        addFieldName_ = args.options()["field"];
        calcType_ = FIELD;
    }
    else if (args.options().found("value"))
    {
        addValueStr_ = args.options()["value"];
        calcType_ = VALUE;
    }
    else
    {
        FatalErrorIn("calcTypes::add::preCalc")
            << "add requires either -field or -value option"
            << nl << exit(FatalError);
    }

    if (args.options().found("resultName"))
    {
        resultName_ = args.options()["resultName"];
    }
}


void Foam::calcTypes::add::calc
(
    const argList& args,
    const Time& runTime,
    const fvMesh& mesh
)
{
    IOobject baseFieldHeader
    (
        baseFieldName_,
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ
    );

    if (baseFieldHeader.headerOk())
    {
        switch (calcType_)
        {
            case FIELD:
            {
                writeAddFields(runTime, mesh, baseFieldHeader);
                break;
            }
            case VALUE:
            {
                writeAddValues(runTime, mesh, baseFieldHeader);
                break;
            }
            default:
            {
                FatalErrorIn("calcTypes::add::calc")
                    << "unknown calcType " << calcType_ << nl
                    << abort(FatalError);
            }
        }
    }
    else
    {
        FatalErrorIn("calcTypes::add::calc")
            << "Unable to read base field: " << baseFieldName_
            << nl << exit(FatalError);
    }
}


// ************************************************************************* //

