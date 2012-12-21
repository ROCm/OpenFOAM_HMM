/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*----------------------------------------------------------------------------*/

#include "MRFSource.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "MRFZone.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(MRFSource, 0);
    addToRunTimeSelectionTable
    (
        basicSource,
        MRFSource,
        dictionary
    );
}

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::MRFSource::initialise()
{
    if (selectionMode_ != smCellZone)
    {
        FatalErrorIn("void Foam::MRFSource::initialise()")
            << "The porosity region must be specified as a cellZone.  Current "
            << "selection mode is " << selectionModeTypeNames_[selectionMode_]
            << exit(FatalError);
    }

    mrfPtr_.reset
    (
        new MRFZone
        (
            name_,
            mesh_,
            coeffs_,
            cellSetName_
        )
    );

    const volVectorField& U = mesh_.lookupObject<volVectorField>(UName_);

    mrfPtr_->correctBoundaryVelocity(const_cast<volVectorField&>(U));

    fieldNames_.setSize(1, UName_);
    applied_.setSize(1, false);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MRFSource::MRFSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    basicSource(name, modelType, dict, mesh),
    mrfPtr_(NULL),
    UName_(coeffs_.lookupOrDefault<word>("UName", "U")),
    rhoName_(coeffs_.lookupOrDefault<word>("rhoName", "rho"))
{
    initialise();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::MRFSource::addSup
(
    fvMatrix<vector>& eqn,
    const label fieldI
)
{
    if (eqn.dimensions() == dimForce)
    {
        const volScalarField& rho =
            mesh_.lookupObject<volScalarField>(rhoName_);

        // use 'true' flag to add to rhs of equation
        mrfPtr_->addCoriolis(rho, eqn, true);
    }
    else
    {
        // use 'true' flag to add to rhs of equation
        mrfPtr_->addCoriolis(eqn, true);
    }
}


void Foam::MRFSource::relativeFlux(surfaceScalarField& phi) const
{
    mrfPtr_->relativeFlux(phi);
}


void Foam::MRFSource::relativeFlux
(
    const surfaceScalarField& rho,
    surfaceScalarField& phi
) const
{
    mrfPtr_->relativeFlux(rho, phi);
}


void Foam::MRFSource::absoluteFlux(surfaceScalarField& phi) const
{
    mrfPtr_->absoluteFlux(phi);
}


void Foam::MRFSource::absoluteFlux
(
    const surfaceScalarField& rho,
    surfaceScalarField& phi
) const
{
    mrfPtr_->absoluteFlux(rho, phi);
}


void Foam::MRFSource::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);
}


bool Foam::MRFSource::read(const dictionary& dict)
{
    if (basicSource::read(dict))
    {
        coeffs_.readIfPresent("UName", UName_);
        coeffs_.readIfPresent("rhoName", rhoName_);

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
