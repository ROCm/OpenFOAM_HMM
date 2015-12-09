/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd.
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

#include "valueAverage.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(valueAverage, 0);
}


// * * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * //

void Foam::valueAverage::writeFileHeader(Ostream& os) const
{
    writeHeader(os, "Value averages");
    writeCommented(os, "Time");
    forAll(fieldNames_, fieldI)
    {
        writeTabbed(os, fieldNames_[fieldI]);
    }
    os  << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::valueAverage::valueAverage
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
    functionObjectName_(dict.lookup("functionObjectName")),
    fieldNames_(dict.lookup("fields")),
    window_(dict.lookupOrDefault<scalar>("window", -1)),
    totalTime_(fieldNames_.size(), obr_.time().deltaTValue()),
    resetOnRestart_(false),
    log_(true)
{
    if (resetOnRestart_)
    {
        forAll(fieldNames_, fieldI)
        {
            const word& fieldName = fieldNames_[fieldI];

            if (dict.found(fieldName))
            {
                const dictionary& valueDict = dict.subDict(fieldName);
                totalTime_[fieldI] = readScalar(valueDict.lookup("totalTime"));
            }
        }
    }

    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

Foam::valueAverage::~valueAverage()
{}


// * * * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * //

void Foam::valueAverage::read(const dictionary& dict)
{
    if (active_)
    {
        functionObjectFile::read(dict);

        log_ = dict.lookupOrDefault<Switch>("log", true);
    }
}


void Foam::valueAverage::execute()
{
    if (!active_)
    {
        return;
    }

    scalar dt = obr_.time().deltaTValue();

    if (log_) Info<< type() << ": " << name_ << " averages:" << nl;

    file() << obr_.time().timeName();

    DynamicList<label> unprocessedFields(fieldNames_.size());

    forAll(fieldNames_, fieldI)
    {
        const word& fieldName(fieldNames_[fieldI]);
        const word meanName(fieldName + "Mean");

        scalar Dt = totalTime_[fieldI];
        scalar alpha = (Dt - dt)/Dt;
        scalar beta = dt/Dt;

        if (window_ > 0)
        {
            if (Dt - dt >= window_)
            {
                alpha = (window_ - dt)/window_;
                beta = dt/window_;
            }
        }

        bool processed = false;
        calc<scalar>(fieldName, meanName, alpha, beta, processed);
        calc<vector>(fieldName, meanName, alpha, beta, processed);
        calc<sphericalTensor>(fieldName, meanName, alpha, beta, processed);
        calc<symmTensor>(fieldName, meanName, alpha, beta, processed);
        calc<tensor>(fieldName, meanName, alpha, beta, processed);

        if (!processed)
        {
            unprocessedFields.append(fieldI);

            if (writeToFile())
            {
                file() << tab << "n/a";
            }
        }

        totalTime_[fieldI] += dt;
    }

    file()<< endl;

    if (unprocessedFields.size())
    {
        WarningInFunction
            << "From function object: " << functionObjectName_ << nl
            << "Unprocessed fields:" << nl;

        forAll(unprocessedFields, i)
        {
            label fieldI = unprocessedFields[i];
            Info<< "        " << fieldNames_[fieldI] << nl;
        }
        Info<< endl;
    }

    if (log_) Info<< endl;
}


void Foam::valueAverage::end()
{
    // Do nothing
}


void Foam::valueAverage::timeSet()
{
    // Do nothing
}


void Foam::valueAverage::write()
{
    // Do nothing
}


// ************************************************************************* //
