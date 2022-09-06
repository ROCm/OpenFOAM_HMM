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

#include "equalBinWidth.H"
#include "histogramModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace histogramModels
{
    defineTypeNameAndDebug(equalBinWidth, 0);
    addToRunTimeSelectionTable(histogramModel, equalBinWidth, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::histogramModels::equalBinWidth::equalBinWidth
(
    const word& name,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    histogramModel(name, mesh, dict),
    nBins_(0),
    min_(GREAT),
    max_(-GREAT)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::histogramModels::equalBinWidth::read(const dictionary& dict)
{
    if (!histogramModel::read(dict))
    {
        return false;
    }

    nBins_ = dict.getScalar("nBins");

    if (nBins_ < 1)
    {
        FatalIOErrorInFunction(dict)
            << "Number of histogram bins = " << nBins_
            << " cannot be negative or zero."
            << abort(FatalIOError);
    }

    min_ = dict.getOrDefault<scalar>("min", GREAT);
    max_ = dict.getOrDefault<scalar>("max", -GREAT);

    return true;
}


bool Foam::histogramModels::equalBinWidth::write(const bool log)
{
    // Retrieve operand field
    const volScalarField& field = histogramModel::getOrReadField(fieldName());

    // Determine min and max from the operand field
    // if the user did not provide any min or max
    scalar histMax = max_;
    scalar histMin = min_;

    if (max_ == -GREAT)
    {
        histMax = max(field).value();

        if (min_ == GREAT)
        {
            histMin = min(field).value();
        }
        if (log)
        {
            Info<< "    Determined histogram bounds from field"
                << " min/max(" << fieldName() << ") = "
                << histMin << ' ' << histMax << endl;
        }
    }
    else if (min_ == GREAT)
    {
        histMin = 0;
    }

    if (histMax < histMin)
    {
        FatalErrorInFunction
            << "Histogram minimum = " << histMin
            << ", cannot be larger than histogram maximum = " << histMax
            << exit(FatalError);
    }


    // Calculate the mid-points of bins for the graph axis
    pointField binMidPoints(nBins_, Zero);
    const scalar delta = (histMax - histMin)/nBins_;

    {
        scalar x = histMin + 0.5*delta;
        for (point& p : binMidPoints)
        {
            p.x() = x;
            x += delta;
        }
    }


    // Calculate the histogram data
    scalarField dataNormalised(nBins_, Zero);
    labelField dataCount(nBins_, Zero);
    const scalarField& V = mesh().V();

    forAll(field, celli)
    {
        const label bini = (field[celli] - histMin)/delta;
        if (bini >= 0 && bini < nBins_)
        {
            dataNormalised[bini] += V[celli];
            dataCount[bini]++;
        }
    }
    Pstream::listCombineGather(dataNormalised, plusEqOp<scalar>());
    Pstream::listCombineGather(dataCount, plusEqOp<label>());


    // Write histogram data
    histogramModel::write(dataNormalised, dataCount, mag(binMidPoints));

    return true;
}


// ************************************************************************* //
