/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::histogram::writeGraph
(
    const coordSet& coords,
    const word& fieldName,
    const scalarField& normalizedValues,
    const scalarField& absoluteValues
) const
{
    fileName outputPath = baseTimeDir();
    mkDir(outputPath);
    OFstream graphFile
    (
        outputPath
       /formatterPtr_().getFileName
        (
            coords,
            wordList(1, fieldName)
        )
    );

    Log << "    Writing histogram of " << fieldName
        << " to " << graphFile.name() << endl;

    wordList fieldNames(2);
    fieldNames[0] = fieldName;
    fieldNames[1] = fieldName + "Count";
    List<const scalarField*> yPtrs(2);
    yPtrs[0] = &normalizedValues;
    yPtrs[1] = &absoluteValues;
    formatterPtr_().write(coords, fieldNames, yPtrs, graphFile);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::histogram::histogram
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name),
    max_(-GREAT),
    min_(GREAT)
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

    const word format(dict.get<word>("setFormat"));
    formatterPtr_ = writer<scalar>::New(format);

    return true;
}


bool Foam::functionObjects::histogram::execute()
{
    return true;
}


bool Foam::functionObjects::histogram::write()
{
    Log << type() << " " << name() << " write:" << nl;

    autoPtr<volScalarField> fieldPtr;
    if (obr_.foundObject<volScalarField>(fieldName_))
    {
        Log << "    Looking up field " << fieldName_ << endl;
    }
    else
    {
        Log << "    Reading field " << fieldName_ << endl;
        fieldPtr.reset
        (
            new volScalarField
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
            )
        );
    }

    const volScalarField& field =
    (
         fieldPtr.valid()
       ? fieldPtr()
       : obr_.lookupObject<volScalarField>(fieldName_)
    );


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
    pointField xBin(nBins_);
    const scalar delta = (histMax- histMin)/nBins_;

    scalar x = histMin + 0.5*delta;
    forAll(xBin, i)
    {
        xBin[i] = point(x, 0, 0);
        x += delta;
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

            const coordSet coords
            (
                fieldName_,
                "x",
                xBin,
                mag(xBin)
            );


            // Convert count field from labelField to scalarField
            scalarField count(dataCount.size());
            forAll(count, i)
            {
                count[i] = 1.0*dataCount[i];
            }

            writeGraph(coords, fieldName_, dataNormalized, count);
        }
    }

    return true;
}


// ************************************************************************* //
