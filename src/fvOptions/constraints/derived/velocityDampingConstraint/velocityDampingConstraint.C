/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
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

#include "velocityDampingConstraint.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMatrix.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(velocityDampingConstraint, 0);
    addToRunTimeSelectionTable
    (
        option,
        velocityDampingConstraint,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::velocityDampingConstraint::addDamping(fvMatrix<vector>& eqn)
{
    // Note: we want to add
    //      deltaU/deltaT
    // where deltaT is a local time scale:
    //      U/(cbrt of volume)
    // Since directly manipulating the diagonal we multiply by volume.

    const scalarField& vol = mesh_.V();
    const volVectorField& U = eqn.psi();
    scalarField& diag = eqn.diag();

    label nDamped = 0;

    forAll(U, cellI)
    {
        const scalar magU = mag(U[cellI]);
        if (magU > UMax_)
        {
            const scalar scale = sqr(Foam::cbrt(vol[cellI]));

            diag[cellI] += scale*(magU-UMax_);

            ++nDamped;
        }
    }

    reduce(nDamped, sumOp<label>());

    Info<< type() << " " << name_ << " damped "
        << nDamped << " ("
        << 100*scalar(nDamped)/mesh_.globalData().nTotalCells()
        << "%) of cells" << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::velocityDampingConstraint::velocityDampingConstraint
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::cellSetOption(name, modelType, dict, mesh)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::velocityDampingConstraint::constrain
(
    fvMatrix<vector>& eqn,
    const label fieldi
)
{
    addDamping(eqn);
}


void Foam::fv::velocityDampingConstraint::writeData(Ostream& os) const
{
    dict_.writeEntry(name_, os);
}


bool Foam::fv::velocityDampingConstraint::read(const dictionary& dict)
{
    if (fv::cellSetOption::read(dict))
    {
        coeffs_.readEntry("UMax", UMax_);

        if (!coeffs_.readIfPresent("UNames", fieldNames_))
        {
            fieldNames_.resize(1);
            fieldNames_.first() = coeffs_.getOrDefault<word>("U", "U");
        }

        fv::option::resetApplied();

        return true;
    }

    return false;
}


// ************************************************************************* //
