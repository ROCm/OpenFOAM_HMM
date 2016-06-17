/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd
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

#include "reactionsSensitivityAnalysis.H"
#include "dictionary.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class chemistryType>
void Foam::reactionsSensitivityAnalysis<chemistryType>::createFileNames()
{
    if (writeToFile() && !prodFilePtr_.valid())
    {
        prodFilePtr_ = createFile("production");
        writeHeader(prodFilePtr_(), "production");
        writeFileHeader(prodFilePtr_());

        consFilePtr_ = createFile("consumption");
        writeHeader(consFilePtr_(), "consumption");
        writeFileHeader(consFilePtr_());

        prodIntFilePtr_ = createFile("productionInt");
        writeHeader(prodIntFilePtr_(), "productionInt");
        writeFileHeader(prodIntFilePtr_());

        consIntFilePtr_ = createFile("consumptionInt");
        writeHeader(consIntFilePtr_(), "consumptionInt");
        writeFileHeader(consIntFilePtr_());
    }
}


template<class chemistryType>
void Foam::reactionsSensitivityAnalysis<chemistryType>::writeFileHeader
(
    OFstream& os
)
{
    writeCommented(os, "Reaction");

    forAll(speciesNames_, k)
    {
        os << tab << speciesNames_[k] << tab;
    }

    os << nl << endl;
}


template<class chemistryType>
void Foam::reactionsSensitivityAnalysis<chemistryType>::calculateSpeciesRR
(
    const basicChemistryModel& basicChemistry
)
{

    tmp<DimensionedField<scalar, volMesh>> RRt
    (
        new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "RR",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimMass/dimVolume/dimTime, 0.0)
        )
    );

    DimensionedField<scalar, volMesh>& RR = RRt.ref();

    scalar dt = mesh_.time().deltaT().value();

    endTime_ += dt;

    forAll(production_, specieI)
    {
        forAll(production_[specieI], reactionI)
        {
            RR = basicChemistry.calculateRR(reactionI, specieI);

            if (RR[0] > 0.0)
            {
                production_[specieI][reactionI] = RR[0];
                productionInt_[specieI][reactionI] =+ dt*RR[0];
            }
            else if (RR[0] < 0.0)
            {
                consumption_[specieI][reactionI] = RR[0];
                consumptionInt_[specieI][reactionI] =+ dt*RR[0];
            }
            else
            {
                production_[specieI][reactionI] = 0.0;
                consumption_[specieI][reactionI] = 0.0;
            }
        }
    }
}


template<class chemistryType>
void Foam::reactionsSensitivityAnalysis<chemistryType>::writeSpeciesRR()
{

    consFilePtr_() << "time : " << mesh_.time().value() << tab << nl;
    consFilePtr_() << "delta T : "<< mesh_.time().deltaT().value() << nl << nl;
    prodFilePtr_() << "time : " << mesh_.time().value() << tab << nl;
    prodFilePtr_() << "delta T : "<< mesh_.time().deltaT().value() << nl << nl;

    consIntFilePtr_() << "start time : " << startTime_ << tab
            << "end time :" <<  endTime_ << nl;

    prodIntFilePtr_() << "start time : " << startTime_ << tab
            << "end time :" <<  endTime_ << nl;

    for
    (
        label reactionI = 0; reactionI < nReactions_; reactionI++
    )
    {
        consFilePtr_() << reactionI << tab;
        consIntFilePtr_() << reactionI << tab;
        prodFilePtr_() << reactionI << tab;
        prodIntFilePtr_() << reactionI << tab;

        forAll(speciesNames_, i)
        {
            prodFilePtr_() << production_[i][reactionI] << tab;
            consFilePtr_() << consumption_[i][reactionI] << tab;
            prodIntFilePtr_() << productionInt_[i][reactionI] << tab;
            consIntFilePtr_() << consumptionInt_[i][reactionI] << tab;
            consumptionInt_[i][reactionI] = 0.0;
            productionInt_[i][reactionI] = 0.0;
        }
        consFilePtr_() << nl;
        consIntFilePtr_() << nl;
        prodFilePtr_() << nl;
        prodIntFilePtr_() << nl;
    }
    consFilePtr_() << nl << nl;
    consIntFilePtr_() << nl << nl;
    prodFilePtr_() << nl << nl;
    prodIntFilePtr_() << nl << nl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class chemistryType>
Foam::reactionsSensitivityAnalysis<chemistryType>::reactionsSensitivityAnalysis
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    functionObjectFile(obr, name),
    name_(name),
    mesh_(refCast<const fvMesh>(obr)),
    active_(true),
    production_(0),
    consumption_(0),
    productionInt_(0),
    consumptionInt_(0),
    startTime_(0),
    endTime_(0),
    speciesNames_(),
    nReactions_(0),
    prodFilePtr_(),
    consFilePtr_(),
    prodIntFilePtr_(),
    consIntFilePtr_()
{
    read(dict);
    if (mesh_.nCells() != 1)
    {
        FatalErrorInFunction
            << "Function object only applicable to single cell cases "
            << abort(FatalError);
    }

    if (mesh_.foundObject<basicChemistryModel>("chemistryProperties"))
    {
        const chemistryType& chemistry = refCast<const chemistryType>
        (
            mesh_.lookupObject<basicChemistryModel>("chemistryProperties")
        );


        speciesNames_.setSize
        (
            chemistry.thermo().composition().species().size()
        );

        forAll(speciesNames_, i)
        {
            speciesNames_[i] = chemistry.thermo().composition().species()[i];
        }

        nReactions_ = chemistry.nReaction();

        if (production_.size() == 0)
        {
            production_.setSize(speciesNames_.size());
            consumption_.setSize(production_.size());
            productionInt_.setSize(production_.size());
            consumptionInt_.setSize(production_.size());

            forAll(production_, i)
            {
                production_[i].setSize(nReactions_, 0.0);
                consumption_[i].setSize(nReactions_, 0.0);
                productionInt_[i].setSize(nReactions_, 0.0);
                consumptionInt_[i].setSize(nReactions_, 0.0);
            }
        }
    }
    else
    {
         FatalErrorInFunction
            << " Not chemistry model found. "
            << " Object available are : " << mesh_.names()
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class chemistryType>
Foam::reactionsSensitivityAnalysis<chemistryType>::
~reactionsSensitivityAnalysis()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class chemistryType>
void Foam::reactionsSensitivityAnalysis<chemistryType>::read
(
    const dictionary& dict
)
{}


template<class chemistryType>
void Foam::reactionsSensitivityAnalysis<chemistryType>::execute()
{
    createFileNames();

    const basicChemistryModel& chemistry =
        mesh_.lookupObject<basicChemistryModel>
        (
            "chemistryProperties"
        );
    calculateSpeciesRR(chemistry);
}


template<class chemistryType>
void Foam::reactionsSensitivityAnalysis<chemistryType>::end()
{
    // Do nothing - only valid on write
}


template<class chemistryType>
void Foam::reactionsSensitivityAnalysis<chemistryType>::timeSet()
{
    // Do nothing
}


template<class chemistryType>
void Foam::reactionsSensitivityAnalysis<chemistryType>::write()
{
    if (!active_)
    {
        return;
    }

    if (Pstream::master())
    {
        //functionObjectFile::write();

        writeSpeciesRR();

        startTime_ = endTime_;
    }
}


// ************************************************************************* //
