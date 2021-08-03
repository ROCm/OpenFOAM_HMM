/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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
#include "areaFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fa
{
    defineTypeNameAndDebug(limitVelocity, 0);
    addToRunTimeSelectionTable
    (
        option,
        limitVelocity,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fa::limitVelocity::limitVelocity
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvPatch& patch
)
:
    faceSetOption(name, modelType, dict, patch),
    UName_(coeffs_.getOrDefault<word>("U", "U")),
    max_(coeffs_.get<scalar>("max"))
{
    fieldNames_.setSize(1, UName_);
    applied_.setSize(1, false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fa::limitVelocity::read(const dictionary& dict)
{
    if (faceSetOption::read(dict))
    {
        coeffs_.readEntry("max", max_);

        return true;
    }

    return false;
}


void Foam::fa::limitVelocity::correct(areaVectorField& U)
{
    const scalar maxSqrU = sqr(max_);

    vectorField& Uif = U.primitiveFieldRef();

    for (const label facei : faces_)
    {
        const scalar magSqrUi = magSqr(Uif[facei]);

        if (magSqrUi > maxSqrU)
        {
            Uif[facei] *= sqrt(maxSqrU/max(magSqrUi, SMALL));
        }
    }

    // Handle boundaries in the case of 'all'
    if (!faceSetOption::useSubMesh())
    {
        for (faPatchVectorField& Up : U.boundaryFieldRef())
        {
            if (!Up.fixesValue())
            {
                forAll(Up, facei)
                {
                    const scalar magSqrUi = magSqr(Up[facei]);

                    if (magSqrUi > maxSqrU)
                    {
                        Up[facei] *= sqrt(maxSqrU/max(magSqrUi, SMALL));
                    }
                }
            }
        }
    }

    // We've changed internal values so give boundary conditions opportunity
    // to correct.
    U.correctBoundaryConditions();
}


// ************************************************************************* //
