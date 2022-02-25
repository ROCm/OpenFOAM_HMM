/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

#include "histogram.H"
#include "volFields.H"
#include "ListOps.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(histogram, 0);
    addToRunTimeSelectionTable(functionObject, histogram, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::histogram::histogram
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObjects::fvMeshFunctionObject(name, runTime, dict),
    functionObjects::writeFile(obr_, name),
    max_(-GREAT),
    min_(GREAT),
    setWriterPtr_(nullptr)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::histogram::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    dict.readEntry("field", fieldName_);

    max_ = dict.getOrDefault<scalar>("max", -GREAT);
    min_ = dict.getOrDefault<scalar>("min", GREAT);
    dict.readEntry("nBins", nBins_);

    if (nBins_ < 1)
    {
        FatalErrorInFunction
            << "Number of histogram bins = " << nBins_
            << " cannot be negative or zero."
            << abort(FatalError);
    }

    const word writeType(dict.get<word>("setFormat"));

    setWriterPtr_ = coordSetWriter::New
    (
        writeType,
        dict.subOrEmptyDict("formatOptions").optionalSubDict(writeType)
    );

    return true;
}


bool Foam::functionObjects::histogram::execute()
{
    return true;
}


bool Foam::functionObjects::histogram::write()
{
    Log << type() << " " << name() << " write:" << nl;

    tmp<volScalarField> tfield;
    tfield.cref(obr_.cfindObject<volScalarField>(fieldName_));

    if (tfield)
    {
        Log << "    Looking up field " << fieldName_ << endl;
    }
    else
    {
        Log << "    Reading field " << fieldName_ << endl;
        tfield = tmp<volScalarField>::New
        (
            IOobject
            (
                fieldName_,
                mesh_.time().timeName(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh_
        );
    }
    const auto& field = tfield();

    scalar histMax = max_;
    scalar histMin = min_;

    if (max_ == -GREAT)
    {
        // Determine current min and max
        histMax = max(field).value();

        if (min_ == GREAT)
        {
            histMin = min(field).value();
        }
        Log << "    Determined histogram bounds from field"
            << " min/max(" << fieldName_ << ") = "
            << histMin << ' ' << histMax << endl;
    }
    else if (min_ == GREAT)
    {
        histMin = 0;
    }

    // Calculate the mid-points of bins for the graph axis
    pointField xBin(nBins_, Zero);
    const scalar delta = (histMax - histMin)/nBins_;

    {
        scalar x = histMin + 0.5*delta;
        for (point& p : xBin)
        {
            p.x() = x;
            x += delta;
        }
    }

    scalarField dataNormalized(nBins_, Zero);
    labelField dataCount(nBins_, Zero);
    const scalarField& V = mesh_.V();

    forAll(field, celli)
    {
        const label bini = (field[celli] - histMin)/delta;
        if (bini >= 0 && bini < nBins_)
        {
            dataNormalized[bini] += V[celli];
            dataCount[bini]++;
        }
    }

    Pstream::listCombineGather(dataNormalized, plusEqOp<scalar>());
    Pstream::listCombineGather(dataCount, plusEqOp<label>());

    if (Pstream::master())
    {
        const scalar sumData = sum(dataNormalized);

        if (sumData > SMALL)
        {
            dataNormalized /= sumData;

            const coordSet coords(fieldName_, "x", xBin, mag(xBin));

            auto& writer = *setWriterPtr_;

            writer.open
            (
                coords,
                (
                    writeFile::baseTimeDir()
                  / (coords.name() + coordSetWriter::suffix(fieldName_))
                )
            );

            Log << "    Writing histogram of " << fieldName_
                << " to " << writer.path() << endl;

            writer.nFields(2);
            writer.write(fieldName_, dataNormalized);
            writer.write(fieldName_ + "Count", dataCount);

            writer.close(true);
        }
    }

    return true;
}


// ************************************************************************* //
