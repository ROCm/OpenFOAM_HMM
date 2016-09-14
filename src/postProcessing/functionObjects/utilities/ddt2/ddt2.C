/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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

#include "ddt2.H"

#include "volFields.H"
#include "dictionary.H"
#include "FieldFunctions.H"
#include "fvcDdt.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ddt2, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


template<class FieldType>
bool Foam::ddt2::calculate
(
    const fvMesh& mesh,
    bool& done
)
{
    if (done)
    {
        return true; // already done - skip
    }

    done = mesh.foundObject<FieldType>(fieldName_);
    if (!done)
    {
        return false;
    }

    const FieldType& input =
        mesh.lookupObject<FieldType>(fieldName_);

    if (!mesh.foundObject<volScalarField>(resultName_))
    {
        const dimensionSet dims
        (
            mag_
          ? mag(input.dimensions()/dimTime)
          : magSqr(input.dimensions()/dimTime)
        );

        mesh.objectRegistry::store
        (
            new volScalarField
            (
                IOobject
                (
                    resultName_,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE // OR IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    "zero",
                    dims,
                    Zero
                ),
                emptyPolyPatch::typeName
            )
        );
    }

    volScalarField& field = const_cast<volScalarField&>
    (
        mesh.lookupObject<volScalarField>(resultName_)
    );

    if (mag_)
    {
        field = mag(fvc::ddt(input));
    }
    else
    {
        field = magSqr(fvc::ddt(input));
    }

    return done;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ddt2::ddt2
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
    fieldName_("undefined-fieldName"),
    resultName_(word::null),
    log_(true),
    mag_(dict.lookupOrDefault<Switch>("mag", false))
{
    // Check if the available mesh is an fvMesh, otherwise deactivate
    if (!isA<fvMesh>(obr_))
    {
        active_ = false;
        WarningInFunction
            << "No fvMesh available, deactivating." << nl
            << endl;
    }

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::ddt2::~ddt2()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ddt2::read(const dictionary& dict)
{
    if (active_)
    {
        log_.readIfPresent("log", dict);

        dict.lookup("fieldName") >> fieldName_;
        dict.readIfPresent("resultName", resultName_);

        if (resultName_.empty())
        {
            resultName_ =
            (
                word(mag_ ? "mag" : "magSqr")
              + "(ddt(" + fieldName_ + "))"
            );
        }
    }
}


void Foam::ddt2::execute()
{
    if (active_)
    {
        const fvMesh& mesh = refCast<const fvMesh>(obr_);
        bool done = false;

        calculate<volScalarField>(mesh, done);
        calculate<volVectorField>(mesh, done);

        if (!done)
        {
            WarningInFunction
                << "Unprocessed field " << fieldName_ << endl;
        }
    }
}


void Foam::ddt2::end()
{
    // Do nothing
}


void Foam::ddt2::timeSet()
{
    // Do nothing
}


void Foam::ddt2::write()
{
    if (active_)
    {
        if (obr_.foundObject<regIOobject>(resultName_))
        {
            const regIOobject& io =
                obr_.lookupObject<regIOobject>(resultName_);

            if (log_)
            {
                const volScalarField& field = dynamicCast<const volScalarField&>
                (
                    io
                );

                // could add additional statistics here
                Info<< type() << " " << name_
                    << " output: writing field " << field.name()
                    << " average: " << gAverage(field) << endl;
            }

            // could also add option to suppress writing?
            io.write();
        }
    }
}


// ************************************************************************* //
