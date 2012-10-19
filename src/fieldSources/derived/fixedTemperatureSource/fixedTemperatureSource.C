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

#include "fixedTemperatureSource.H"
#include "fvMesh.H"
#include "fvMatrices.H"
#include "basicThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fixedTemperatureSource, 0);
    addToRunTimeSelectionTable
    (
        basicSource,
        fixedTemperatureSource,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedTemperatureSource::fixedTemperatureSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    basicSource(name, modelType, dict, mesh),
    T_(readScalar(coeffs_.lookup("temperature")))
{
    fieldNames_.setSize(1, "energy");
    applied_.setSize(1, false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fixedTemperatureSource::alwaysApply() const
{
    return true;
}


void Foam::fixedTemperatureSource::setValue
(
    fvMatrix<scalar>& eqn,
    const label
)
{
    const basicThermo& thermo =
        mesh_.lookupObject<basicThermo>("thermophysicalProperties");

    if (eqn.psi().name() == thermo.he().name())
    {
        const scalarField Tfield(cells_.size(), T_);
        eqn.setValues(cells_, thermo.he(thermo.p(), Tfield, cells_));
    }
}


void Foam::fixedTemperatureSource::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);
}


bool Foam::fixedTemperatureSource::read(const dictionary& dict)
{
    if (basicSource::read(dict))
    {
        coeffs_.readIfPresent("T", T_);

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
