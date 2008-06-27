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

#include "fieldAverage.H"
#include "volFields.H"
#include "dictionary.H"
#include "Time.H"
#include "IFstream.H"
#include "OFstream.H"

#include "fieldAverageItem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fieldAverage, 0);
}

const Foam::word Foam::fieldAverage::EXT_MEAN = "Mean";
const Foam::word Foam::fieldAverage::EXT_PRIME2MEAN = "Prime2Mean";


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fieldAverage::resetLists(const label nItems)
{
    meanScalarFields_.clear();
    meanScalarFields_.setSize(nItems);

    meanVectorFields_.clear();
    meanVectorFields_.setSize(nItems);

    meanSphericalTensorFields_.clear();
    meanSphericalTensorFields_.setSize(nItems);

    meanSymmTensorFields_.clear();
    meanSymmTensorFields_.setSize(nItems);

    meanTensorFields_.clear();
    meanTensorFields_.setSize(nItems);

    prime2MeanScalarFields_.clear();
    prime2MeanScalarFields_.setSize(nItems);

    prime2MeanSymmTensorFields_.clear();
    prime2MeanSymmTensorFields_.setSize(nItems);

    totalIter_.clear();
    totalIter_.setSize(nItems, 1);

    totalTime_.clear();
    totalTime_.setSize(nItems, obr_.time().deltaT().value());
}


void Foam::fieldAverage::initialise()
{
    // Add mean fields to the field lists
    forAll(faItems_, i)
    {
        const word fieldName = faItems_[i].fieldName();
        if (obr_.foundObject<volScalarField>(fieldName))
        {
            addMeanField<scalar>(i, meanScalarFields_);
        }
        else if (obr_.foundObject<volVectorField>(fieldName))
        {
            addMeanField<vector>(i, meanVectorFields_);
        }
        else if (obr_.foundObject<volSphericalTensorField>(fieldName))
        {
            addMeanField<sphericalTensor>(i, meanSphericalTensorFields_);
        }
        else if (obr_.foundObject<volSymmTensorField>(fieldName))
        {
            addMeanField<symmTensor>(i, meanSymmTensorFields_);
        }
        else if (obr_.foundObject<volTensorField>(fieldName))
        {
            addMeanField<tensor>(i, meanTensorFields_);
        }
        else
        {
            FatalErrorIn("Foam::fieldAverage::initialise()")
                << "Requested field " << faItems_[i].fieldName()
                << " does not exist in the database" << nl
                << exit(FatalError);
        }
    }

    // Add prime-squared mean fields to the field lists
    forAll(faItems_, i)
    {
        if (faItems_[i].prime2Mean())
        {
            const word fieldName = faItems_[i].fieldName();
            if (!faItems_[i].mean())
            {
                FatalErrorIn("Foam::fieldAverage::initialise()")
                    << "To calculate the prime-squared average, the "
                    << "mean average must also be selected for field "
                    << fieldName << nl << exit(FatalError);
            }

            if (obr_.foundObject<volScalarField>(fieldName))
            {
                addPrime2MeanField<scalar>
                (
                    i,
                    meanScalarFields_,
                    prime2MeanScalarFields_
                );
            }
            else if (obr_.foundObject<volVectorField>(fieldName))
            {
                addPrime2MeanField<vector>
                (
                    i,
                    meanVectorFields_,
                    prime2MeanSymmTensorFields_
                );
            }
            else
            {
                FatalErrorIn("Foam::fieldAverage::initialise()")
                    << "prime2Mean average can only be applied to "
                    << "volScalarFields and volVectorFields"
                    << nl << "    Field: " << fieldName << nl
                    << exit(FatalError);
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fieldAverage::fieldAverage
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
    faItems_(dict.lookup("fields")),
    meanScalarFields_(faItems_.size()),
    meanVectorFields_(faItems_.size()),
    meanSphericalTensorFields_(faItems_.size()),
    meanSymmTensorFields_(faItems_.size()),
    meanTensorFields_(faItems_.size()),
    prime2MeanScalarFields_(faItems_.size()),
    prime2MeanSymmTensorFields_(faItems_.size()),
    totalIter_(faItems_.size(), 1),
    totalTime_(faItems_.size(), obr_.time().deltaT().value())
{
    // Check if the available mesh is an fvMesh otherise deactivate
    if (!isA<fvMesh>(obr_))
    {
        active_ = false;
        WarningIn
        (
            "fieldAverage::fieldAverage"
            "("
                "const word&,"
                "const objectRegistry&,"
                "const dictionary&,"
                "const bool"
            ")"
        )   << "No fvMesh available, deactivating."
            << nl << endl;
    }

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fieldAverage::~fieldAverage()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fieldAverage::read(const dictionary& dict)
{
    if (active_)
    {
        faItems_.clear();
        faItems_ = List<fieldAverageItem>(dict.lookup("fields"));

        resetLists(faItems_.size());

        initialise();

        readAveragingProperties();
    }
}


void Foam::fieldAverage::write()
{
    if (active_)
    {
        calcAverages();

        if (obr_.time().outputTime())
        {
            writeAverages();

            writeAveragingProperties();
        }
    }
}


void Foam::fieldAverage::calcAverages()
{
    Info<< "Calculating averages" << nl << endl;
    forAll(faItems_, i)
    {
        totalIter_[i]++;
        totalTime_[i] += obr_.time().deltaT().value();
    }

    addMeanSqrToPrime2Mean<scalar>
    (
        meanScalarFields_,
        prime2MeanScalarFields_
    );
    addMeanSqrToPrime2Mean<vector>
    (
        meanVectorFields_,
        prime2MeanSymmTensorFields_
    );

    calculateMeanFields<scalar>(meanScalarFields_);
    calculateMeanFields<vector>(meanVectorFields_);
    calculateMeanFields<sphericalTensor>(meanSphericalTensorFields_);
    calculateMeanFields<symmTensor>(meanSymmTensorFields_);
    calculateMeanFields<tensor>(meanTensorFields_);

    calculatePrime2MeanFields<scalar>
    (
        meanScalarFields_,
        prime2MeanScalarFields_
    );
    calculatePrime2MeanFields<vector>
    (
        meanVectorFields_,
        prime2MeanSymmTensorFields_
    );
}


void Foam::fieldAverage::writeAverages() const
{
    writeFieldList<scalar>(meanScalarFields_);
    writeFieldList<vector>(meanVectorFields_);
    writeFieldList<sphericalTensor>(meanSphericalTensorFields_);
    writeFieldList<symmTensor>(meanSymmTensorFields_);
    writeFieldList<tensor>(meanTensorFields_);

    writeFieldList<scalar>(prime2MeanScalarFields_);
    writeFieldList<symmTensor>(prime2MeanSymmTensorFields_);
}


void Foam::fieldAverage::writeAveragingProperties() const
{
    IOdictionary propsDict
    (
        IOobject
        (
            "fieldAveragingProperties",
            obr_.time().timeName(),
            "uniform",
            obr_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    );

    forAll(faItems_, i)
    {
        const word fieldName = faItems_[i].fieldName();
        propsDict.add(fieldName, dictionary());
        propsDict.subDict(fieldName).add("totalIter", totalIter_[i]);
        propsDict.subDict(fieldName).add("totalTime", totalTime_[i]);
    }

    propsDict.regIOobject::write();
}


void Foam::fieldAverage::readAveragingProperties()
{
    IFstream propsFile
    (
        obr_.time().path()/obr_.time().timeName()
       /"uniform"/"fieldAveragingProperties"
    );

    if (!propsFile.good())
    {
        return;
    }

    dictionary propsDict(dictionary::null, propsFile);

    forAll(faItems_, i)
    {
        const word& fieldName = faItems_[i].fieldName();
        if (propsDict.found(fieldName))
        {
            dictionary fieldDict(propsDict.subDict(fieldName));

            totalIter_[i] = readLabel(fieldDict.lookup("totalIter"));
            totalTime_[i] = readScalar(fieldDict.lookup("totalTime"));
        }
    }
}


void Foam::fieldAverage::updateMesh(const mapPolyMesh&)
{
    // Do nothing
}


void Foam::fieldAverage::movePoints(const pointField&)
{
    // Do nothing
}


// ************************************************************************* //
