/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "RASModel.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class BasicTurbulenceModel>
void Foam::RASModel<BasicTurbulenceModel>::printCoeffs(const word& type)
{
    if (printCoeffs_)
    {
        Info<< coeffDict_.dictName() << coeffDict_ << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
Foam::RASModel<BasicTurbulenceModel>::RASModel
(
    const word& type,
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
:
    BasicTurbulenceModel
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    RASDict_(this->subOrEmptyDict("RAS")),
    turbulence_(RASDict_.get<Switch>("turbulence")),
    printCoeffs_(RASDict_.getOrDefault<Switch>("printCoeffs", false)),
    coeffDict_(RASDict_.optionalSubDict(type + "Coeffs")),

    kMin_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "kMin",
            RASDict_,
            sqr(dimVelocity),
            SMALL
        )
    ),

    epsilonMin_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "epsilonMin",
            RASDict_,
            kMin_.dimensions()/dimTime,
            SMALL
        )
    ),

    omegaMin_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "omegaMin",
            RASDict_,
            dimless/dimTime,
            SMALL
        )
    )
{
    // Force the construction of the mesh deltaCoeffs which may be needed
    // for the construction of the derived models and BCs
    this->mesh_.deltaCoeffs();
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
Foam::autoPtr<Foam::RASModel<BasicTurbulenceModel>>
Foam::RASModel<BasicTurbulenceModel>::New
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
{
    const IOdictionary modelDict
    (
        IOobject
        (
            IOobject::groupName(propertiesName, alphaRhoPhi.group()),
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false // Do not register
        )
    );

    const dictionary& dict = modelDict.subDict("RAS");

    const word modelType(dict.get<word>("RASModel"));

    Info<< "Selecting RAS turbulence model " << modelType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.found())
    {
        FatalIOErrorInLookup
        (
            dict,
            "RASModel",
            modelType,
            *dictionaryConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<RASModel>
    (
        cstrIter()(alpha, rho, U, alphaRhoPhi, phi, transport, propertiesName)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool Foam::RASModel<BasicTurbulenceModel>::read()
{
    if (BasicTurbulenceModel::read())
    {
        RASDict_ <<= this->subDict("RAS");
        RASDict_.readEntry("turbulence", turbulence_);

        coeffDict_ <<= RASDict_.optionalSubDict(type() + "Coeffs");

        kMin_.readIfPresent(RASDict_);
        epsilonMin_.readIfPresent(RASDict_);
        omegaMin_.readIfPresent(RASDict_);

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
void Foam::RASModel<BasicTurbulenceModel>::correct()
{
    BasicTurbulenceModel::correct();
}


// ************************************************************************* //
