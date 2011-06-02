/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2011 OpenCFD Ltd.
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


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

noPyrolysis::noPyrolysis(const word& modelType, const fvMesh& mesh)
:
    pyrolysisModel(mesh)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

noPyrolysis::~noPyrolysis()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const volScalarField& noPyrolysis::rho() const
{
    FatalErrorIn("const volScalarField& noPyrolysis::rho() const")
        << "rho field not available for " << type() << abort(FatalError);
    return volScalarField::null();
}


const volScalarField& noPyrolysis::T() const
{
    FatalErrorIn("const volScalarField& noPyrolysis::T() const")
        << "T field not available for " << type() << abort(FatalError);
    return volScalarField::null();
}


const tmp<volScalarField> noPyrolysis::Cp() const
{
    FatalErrorIn("const tmp<volScalarField>& noPyrolysis::Cp() const")
        << "Cp field not available for " << type() << abort(FatalError);

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "noPyrolysis::Cp",
                time().timeName(),
                primaryMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            primaryMesh(),
            dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0)
        )
    );
}


const volScalarField& noPyrolysis::kappa() const
{
    FatalErrorIn("const volScalarField& noPyrolysis::kappa() const")
        << "kappa field not available for " << type() << abort(FatalError);
    return volScalarField::null();
}


const volScalarField& noPyrolysis::K() const
{
    FatalErrorIn("const volScalarField& noPyrolysis::K() const")
        << "K field not available for " << type() << abort(FatalError);
    return volScalarField::null();
}


const surfaceScalarField& noPyrolysis::phiGas() const
{
    FatalErrorIn("const volScalarField& noPyrolysis::phiGas() const")
        << "phiGas field not available for " << type() << abort(FatalError);
    return surfaceScalarField::null();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
