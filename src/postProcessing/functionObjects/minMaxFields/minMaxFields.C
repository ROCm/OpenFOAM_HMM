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

#include "minMaxFields.H"
#include "volFields.H"
#include "dictionary.H"
#include "Time.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(minMaxFields, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::minMaxFields::minMaxFields
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    active_(true),
    log_(false),
    fieldSet_(),
    minMaxFieldsFilePtr_(NULL)
{
    // Check if the available mesh is an fvMesh otherise deactivate
    if (!isA<fvMesh>(obr_))
    {
        active_ = false;
        WarningIn
        (
            "minMaxFields::minMaxFields"
            "(const objectRegistry& obr, const dictionary& dict)"
        )   << "No fvMesh available, deactivating."
            << endl;
    }

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::minMaxFields::~minMaxFields()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::minMaxFields::read(const dictionary& dict)
{
    if (active_)
    {
        log_ = dict.lookupOrDefault<Switch>("log", false);

        dict.lookup("fields") >> fieldSet_;
    }
}


void Foam::minMaxFields::makeFile()
{
    // Create the minMaxFields file if not already created
    if (!minMaxFieldsFilePtr_.valid())
    {
        if (debug)
        {
            Info<< "Creating minMaxFields file." << endl;
        }

        // File update
        if (Pstream::master())
        {
            fileName minMaxFieldsDir;
            if (Pstream::parRun())
            {
                // Put in undecomposed case (Note: gives problems for
                // distributed data running)
                minMaxFieldsDir =
                    obr_.time().path()/".."/name_/obr_.time().timeName();
            }
            else
            {
                minMaxFieldsDir =
                    obr_.time().path()/name_/obr_.time().timeName();
            }

            // Create directory if does not exist.
            mkDir(minMaxFieldsDir);

            // Open new file at start up
            minMaxFieldsFilePtr_.reset
            (
                new OFstream(minMaxFieldsDir/(type() + ".dat"))
            );

            // Add headers to output data
            writeFileHeader();
        }
    }
}


void Foam::minMaxFields::writeFileHeader()
{
    if (minMaxFieldsFilePtr_.valid())
    {
        minMaxFieldsFilePtr_()
            << "# Time" << tab << "field" << tab << "min" << tab << "max"
            << endl;
    }
}


void Foam::minMaxFields::execute()
{
    // Do nothing - only valid on write
}

void Foam::minMaxFields::write()
{
    if (active_)
    {
        // Create the minMaxFields file if not already created
        makeFile();

        forAll(fieldSet_, fieldI)
        {
            calcMinMaxFields<scalar>(fieldSet_[fieldI]);
            calcMinMaxFields<vector>(fieldSet_[fieldI]);
            calcMinMaxFields<sphericalTensor>(fieldSet_[fieldI]);
            calcMinMaxFields<symmTensor>(fieldSet_[fieldI]);
            calcMinMaxFields<tensor>(fieldSet_[fieldI]);
        }
    }
}


template<>
void Foam::minMaxFields::calcMinMaxFields<Foam::scalar>
(
    const word& fieldName
)
{
    if (obr_.foundObject<volScalarField>(fieldName))
    {
        const scalarField& field = obr_.lookupObject<scalarField>(fieldName);
        scalar minValue = min(field);
        scalar maxValue = max(field);

        reduce(minValue, minOp<scalar>());
        reduce(maxValue, maxOp<scalar>());

        if (Pstream::master())
        {
            minMaxFieldsFilePtr_() << obr_.time().value() << tab
                << fieldName << tab << minValue << tab << maxValue << endl;

            if (log_)
            {
                Info<< "minMaxFields output:" << nl
                    << "    min(" << fieldName << ") = " << minValue << nl
                    << "    max(" << fieldName << ") = " << maxValue << nl
                    << endl;
            }
        }
    }
}


// ************************************************************************* //
