/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "noPyrolysis.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace pyrolysisModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(noPyrolysis, 0);
addToRunTimeSelectionTable(pyrolysisModel, noPyrolysis, mesh);
addToRunTimeSelectionTable(pyrolysisModel, noPyrolysis, dictionary);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool noPyrolysis::read()
{
    if (pyrolysisModel::read())
    {
        // no additional info to read
        return true;
    }
    else
    {
        return false;
    }
}


bool noPyrolysis::read(const dictionary& dict)
{
    if (pyrolysisModel::read(dict))
    {
        // no additional info to read
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

noPyrolysis::noPyrolysis(const word& modelType, const fvMesh& mesh)
:
    pyrolysisModel(modelType, mesh),
    solidChemistry_(solidChemistryModel::New(regionMesh())),
    solidThermo_(solidChemistry_->solidThermo())
{}


noPyrolysis::noPyrolysis
(
    const word& modelType,
    const fvMesh& mesh,
    const dictionary& dict
):
    pyrolysisModel(modelType, mesh, dict),
    solidChemistry_(solidChemistryModel::New(regionMesh())),
    solidThermo_(solidChemistry_->solidThermo())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

noPyrolysis::~noPyrolysis()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void noPyrolysis::preEvolveRegion()
{
    //Do nothing
}


void noPyrolysis::evolveRegion()
{
    //Do nothing
}

const volScalarField& noPyrolysis::rho() const
{
    return (solidThermo_.rho());
}


const volScalarField& noPyrolysis::T() const
{
    return (solidThermo_.T());
}


const tmp<volScalarField> noPyrolysis::Cp() const
{
    return (solidThermo_.Cp());
}


const volScalarField& noPyrolysis::kappa() const
{
    return (solidThermo_.kappa());
}


const volScalarField& noPyrolysis::K() const
{
    return (solidThermo_.K());
}


const surfaceScalarField& noPyrolysis::phiGas() const
{
    FatalErrorIn("const volScalarField& noPyrolysis::phiGas() const")
        << "phiGas field not available for " << type() << abort(FatalError);
    return surfaceScalarField::null();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace pyrolysisFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
