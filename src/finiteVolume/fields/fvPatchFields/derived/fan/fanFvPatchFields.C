/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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

#include "fanFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<>
void Foam::fanFvPatchField<Foam::scalar>::calcFanJump()
{
    if (this->cyclicPatch().owner())
    {
        const surfaceScalarField& phi =
            db().lookupObject<surfaceScalarField>(phiName_);

        const fvsPatchField<scalar>& phip =
            patch().patchField<surfaceScalarField, scalar>(phi);

        scalarField Un(max(phip/patch().magSf(), scalar(0)));
        if (uniformJump_)
        {
            scalar area = gSum(patch().magSf());
            Un = gSum(Un*patch().magSf())/area;

            if (nonDimensional_)
            {
                // Create an non-dimensional velocity
                Un =
                    120.0*Un/pow3(constant::mathematical::pi)
                  / dm_/rpm_;
            }
        }

        if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
        {
            Un /= patch().lookupPatchField<volScalarField, scalar>(rhoName_);
        }

        if (nonDimensional_)
        {
            scalarField deltap(this->jumpTable_->value(Un));

            // Convert non-dimensional deltap from curve into deltaP
            scalarField pdFan
            (
                deltap*pow4(constant::mathematical::pi)*sqr(dm_*rpm_)/1800.0
            );

            this->setJump(pdFan);
        }
        else
        {
            this->setJump(jumpTable_->value(Un));
        }

        this->relax();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeTemplatePatchTypeField(scalar, fan);
}


// ************************************************************************* //
