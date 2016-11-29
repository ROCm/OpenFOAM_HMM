/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 OpenCFD Ltd.
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

#include "binned.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace distributionModels
    {
        defineTypeNameAndDebug(binned, 0);
        addToRunTimeSelectionTable(distributionModel, binned, dictionary);
    }
}


const char* Foam::distributionModels::binned::header =
    "minValue maxValue (binBoundaries) (values)";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distributionModels::binned::binned
(
    const dictionary& dict,
    cachedRandom& rndGen
)
:
    distributionModel(typeName, dict, rndGen),
    minValue_(readScalar(distributionModelDict_.lookup("minValue"))),
    maxValue_(readScalar(distributionModelDict_.lookup("maxValue"))),
    meanValue_(0),
    binBoundaries_(distributionModelDict_.lookup("binBoundaries")),
    values_(distributionModelDict_.lookup("values"))
{
    if (minValue_ < 0)
    {
        FatalErrorInFunction
            << "Minimum value must be greater than zero. "
            << "Supplied minValue = " << minValue_
            << exit(FatalError);
    }

    if (maxValue_ < minValue_)
    {
        FatalErrorInFunction
            << "Maximum value is smaller than the minimum value:"
            << "    maxValue = " << maxValue_ << ", minValue = " << minValue_
            << exit(FatalError);
    }

    if (binBoundaries_.size() != values_.size() + 1)
    {
        FatalErrorInFunction
            << "Number of bin boundaries must be number of values + 1"
            << "    number of binBoundaries : " << binBoundaries_.size() << nl
            << "    number of values        : " << values_.size() << nl
            << exit(FatalError);
    }
}


Foam::distributionModels::binned::binned
(
    const scalarField& values,
    const scalar binWidth,
    cachedRandom& rndGen
)
:
    distributionModel(typeName, dictionary::null, rndGen),
    minValue_(0),
    maxValue_(0),
    meanValue_(0),
    binBoundaries_(),
    values_()
{
    label bin0 = floor(min(values)/binWidth);
    label bin1 = ceil(max(values)/binWidth);
    label nBin = bin1 - bin0;

    minValue_ = bin0*binWidth;
    maxValue_ = bin1*binWidth;

    if (nBin == 0)
    {
        WarningInFunction
            << "Data cannot be binned - zero bins generated" << nl
            << "   Bin width : " << binWidth << nl
            << "   Values    : " << values
            << endl;

        return;
    }

    // Create the bins
    binBoundaries_.setSize(nBin + 1);
    forAll(binBoundaries_, binI)
    {
        binBoundaries_[binI] = (bin0 + binI)*binWidth;
    }

    // Bin and normalise the data
    values_.setSize(nBin, 0);
    forAll(values, i)
    {
        label binI = floor(values[i]/binWidth) - bin0;
        values_[binI]++;
    }
    values_ /= sum(values_);

    // Calculate the mean value
    scalar sum = 0;
    label valueI = 0;
    forAll(values_, i)
    {
        sum += values_[i];
        if (sum > 0.5)
        {
            valueI = i;
            break;
        }
    }

    if (valueI == 0)
    {
        meanValue_ = binBoundaries_[valueI];
    }
    else
    {
        meanValue_ =
            binBoundaries_[valueI]
          + (0.5 - binBoundaries_[valueI])
          / (binBoundaries_[valueI + 1] - binBoundaries_[valueI]);
    }
}


Foam::distributionModels::binned::binned(const binned& p)
:
    distributionModel(p),
    minValue_(p.minValue_),
    maxValue_(p.maxValue_),
    binBoundaries_(p.binBoundaries_),
    values_(p.values_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::distributionModels::binned::~binned()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::distributionModels::binned::sample() const
{
    scalar y = rndGen_.sample01<scalar>();

    scalar sum = 0;
    forAll(values_, i)
    {
        sum += values_[i];
        if (sum > y)
        {
            scalar binMin = binBoundaries_[i];
            scalar binMax = binBoundaries_[i + 1];
            return binMin + (1 - (sum - y)/values_[i])*(binMax - binMin);
        }
    }

    return maxValue_;
}


Foam::scalar Foam::distributionModels::binned::minValue() const
{
    return minValue_;
}


Foam::scalar Foam::distributionModels::binned::maxValue() const
{
    return maxValue_;
}


Foam::scalar Foam::distributionModels::binned::meanValue() const
{
    return meanValue_;
}


void Foam::distributionModels::binned::readData(Istream& is)
{
//    distributionModel::readData(is);
    is  >> minValue_
        >> maxValue_
        >> binBoundaries_
        >> values_;
}


void Foam::distributionModels::binned::writeData(Ostream& os) const
{
//    distributionModel::writeData(os);
    os  << minValue_ << token::SPACE
        << maxValue_ << token::SPACE
        << binBoundaries_ << token::SPACE
        << values_;
}


Foam::dictionary Foam::distributionModels::binned::writeDict
(
    const word& dictName
) const
{
//    dictionary dict = distributionModel::writeDict(dictName);
    dictionary dict(dictName);

    dict.add("minValue", minValue_);
    dict.add("maxValue", maxValue_);
    dict.add("binBoundaries", binBoundaries_);
    dict.add("values", values_);

    return dict;
}


void Foam::distributionModels::binned::readDict(const dictionary& dict)
{
//    distributionModel::readDict(dict);

    dict.lookup("minValue") >> minValue_;
    dict.lookup("maxValue") >> maxValue_;
    dict.lookup("binBoundaries") >> binBoundaries_;
    dict.lookup("values") >> values_;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const distributionModels::binned& b)
{
    os.check("Foam::Ostream& Foam::operator<<(Ostream&, const binned&)");

    b.writeData(os);
    return os;
}


Foam::Istream&  Foam::operator>>(Istream& is, distributionModels::binned& b)
{
    is.check("Foam::Istream& Foam::operator>>(Istream&, binned&)");

    b.readData(is);
    return is;
}


// ************************************************************************* //
