/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
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

template<class RhoFieldType>
void Foam::porosityModels::pressureDrop::apply
(
    scalarField& Udiag,
    vectorField& Usource,
    const scalarField& V,
    const RhoFieldType& rho,
    const vectorField& U,
    const scalar rhoScale
) const
{
    // local-to-global transformation tensor
    const tensor& E = coordSys_.R().R();

    forAll(cellZoneIds_, zoneI)
    {
        const labelList& cells = mesh_.cellZones()[cellZoneIds_[zoneI]];

        forAll(cells, i)
        {
            const label cellI = cells[i];
            const scalar magU = mag(U[cellI]);
            const scalar mDot = rho[cellI]*magU;
            const scalar dp = mDotvsDp_->value(mDot);
            const tensor Cd = E*dp/(lRef_*magU*rhoScale + ROOTVSMALL);
            const scalar isoCd = tr(Cd);

            Udiag[cellI] += V[cellI]*isoCd;
            Usource[cellI] -= V[cellI]*((Cd - I*isoCd) & U[cellI]);
        }
    }
}


template<class RhoFieldType>
void Foam::porosityModels::pressureDrop::apply
(
    tensorField& AU,
    const RhoFieldType& rho,
    const vectorField& U,
    const scalar rhoScale
) const
{
    // local-to-global transformation tensor
    const tensor& E = coordSys_.R().R();

    forAll(cellZoneIds_, zoneI)
    {
        const labelList& cells = mesh_.cellZones()[cellZoneIds_[zoneI]];

        forAll(cells, i)
        {
            const label cellI = cells[i];
            const scalar magU = mag(U[cellI]);
            const scalar mDot = rho[cellI]*magU;
            const scalar dp = mDotvsDp_->value(mDot);

            AU[cellI] += E*dp/(lRef_*magU*rhoScale + ROOTVSMALL);
        }
    }
}


// ************************************************************************* //
