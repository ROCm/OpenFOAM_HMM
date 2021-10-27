/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 OpenFOAM Foundation
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

#include "limitVelocity.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(limitVelocity, 0);
    addToRunTimeSelectionTable(option, limitVelocity, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::limitVelocity::limitVelocity
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::cellSetOption(name, modelType, dict, mesh),
    UName_(coeffs_.getOrDefault<word>("U", "U")),
    max_(coeffs_.get<scalar>("max"))
{
    fieldNames_.resize(1, UName_);
    fv::option::resetApplied();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::limitVelocity::read(const dictionary& dict)
{
    if (fv::cellSetOption::read(dict))
    {
        coeffs_.readEntry("max", max_);

        return true;
    }

    return false;
}


void Foam::fv::limitVelocity::correct(volVectorField& U)
{
    const scalar maxSqrU = sqr(max_);

    vectorField& Uif = U.primitiveFieldRef();

    for (const label celli : cells_)
    {
        const scalar magSqrUi = magSqr(Uif[celli]);

        if (magSqrUi > maxSqrU)
        {
            Uif[celli] *= sqrt(maxSqrU/magSqrUi);
        }
    }

    // Handle boundaries in the case of 'all'
    if (!cellSetOption::useSubMesh())
    {
        for (fvPatchVectorField& Up : U.boundaryFieldRef())
        {
            if (!Up.fixesValue())
            {
                forAll(Up, facei)
                {
                    const scalar magSqrUi = magSqr(Up[facei]);

                    if (magSqrUi > maxSqrU)
                    {
                        Up[facei] *= sqrt(maxSqrU/magSqrUi);
                    }
                }
            }
        }
    }

    // We've changed internal values so give
    // boundary conditions opportunity to correct
    U.correctBoundaryConditions();
}


// ************************************************************************* //
