/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

#include "limitTemperature.H"
#include "fvMesh.H"
#include "basicThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(limitTemperature, 0);
    addToRunTimeSelectionTable(option, limitTemperature, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::limitTemperature::limitTemperature
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::cellSetOption(name, modelType, dict, mesh),
    Tmin_(coeffs_.get<scalar>("min")),
    Tmax_(coeffs_.get<scalar>("max")),
    phase_(coeffs_.getOrDefault<word>("phase", word::null))
{
    // Set the field name to that of the energy
    // field from which the temperature is obtained
    const auto& thermo =
        mesh_.lookupObject<basicThermo>
        (
            IOobject::groupName(basicThermo::dictName, phase_)
        );

    fieldNames_.resize(1, thermo.he().name());

    fv::option::resetApplied();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::limitTemperature::read(const dictionary& dict)
{
    if (fv::cellSetOption::read(dict))
    {
        coeffs_.readEntry("min", Tmin_);
        coeffs_.readEntry("max", Tmax_);

        if (Tmax_ < Tmin_)
        {
            FatalIOErrorInFunction(dict)
                << "Minimum temperature limit cannot exceed maximum limit" << nl
                << "min = " << Tmin_ << nl
                << "max = " << Tmax_
                << exit(FatalIOError);
        }

        if (Tmin_ < 0)
        {
            FatalIOErrorInFunction(dict)
                << "Minimum temperature limit cannot be negative" << nl
                << "min = " << Tmin_
                << exit(FatalIOError);
        }

        return true;
    }

    return false;
}


void Foam::fv::limitTemperature::correct(volScalarField& he)
{
    const auto& thermo =
        mesh_.lookupObject<basicThermo>
        (
            IOobject::groupName(basicThermo::dictName, phase_)
        );

    scalarField Tmin(cells_.size(), Tmin_);
    scalarField Tmax(cells_.size(), Tmax_);

    scalarField heMin(thermo.he(thermo.p(), Tmin, cells_));
    scalarField heMax(thermo.he(thermo.p(), Tmax, cells_));

    scalarField& hec = he.primitiveFieldRef();

    const scalarField& T = thermo.T();

    scalar Tmin0 = min(T);
    scalar Tmax0 = max(T);

    label nOverTmax = 0;
    label nLowerTmin = 0;

    forAll(cells_, i)
    {
        const label celli = cells_[i];
        if (hec[celli] < heMin[i])
        {
            nLowerTmin++;
        }
        else if (hec[celli] > heMax[i])
        {
            nOverTmax++;
        }
        hec[celli]= max(min(hec[celli], heMax[i]), heMin[i]);
    }

    reduce(nOverTmax, sumOp<label>());
    reduce(nLowerTmin, sumOp<label>());

    reduce(Tmin0, minOp<scalar>());
    reduce(Tmax0, maxOp<scalar>());

    Info<< type() << " " << name_ << " Lower limited "
        << nLowerTmin << " ("
        << 100*scalar(nLowerTmin)/mesh_.globalData().nTotalCells()
        << "%) of cells" << endl;

    Info<< type() << " " << name_ << " Upper limited "
        << nOverTmax << " ("
        << 100*scalar(nOverTmax)/mesh_.globalData().nTotalCells()
        << "%) of cells" << endl;

    Info<< type() << " " << name_ << " Unlimited Tmax " << Tmax0 << nl
        <<  "Unlimited Tmin " << Tmin0 << endl;

    // Handle boundaries in the case of 'all'
    if (!cellSetOption::useSubMesh())
    {
        volScalarField::Boundary& bf = he.boundaryFieldRef();

        forAll(bf, patchi)
        {
            fvPatchScalarField& hep = bf[patchi];

            if (!hep.fixesValue())
            {
                const scalarField& pp = thermo.p().boundaryField()[patchi];

                scalarField Tminp(pp.size(), Tmin_);
                scalarField Tmaxp(pp.size(), Tmax_);

                scalarField heMinp(thermo.he(pp, Tminp, patchi));
                scalarField heMaxp(thermo.he(pp, Tmaxp, patchi));

                forAll(hep, facei)
                {
                    hep[facei] =
                        max(min(hep[facei], heMaxp[facei]), heMinp[facei]);
                }
            }
        }
    }

    // We've changed internal values so give
    // boundary conditions opportunity to correct
    he.correctBoundaryConditions();
}


// ************************************************************************* //
