/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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

#include "thermoSingleLayer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const Foam::SLGThermo&
Foam::surfaceFilmModels::thermoSingleLayer::thermo() const
{
    return thermo_;
}


inline Foam::label Foam::surfaceFilmModels::thermoSingleLayer::liquidId() const
{
    return liquidId_;
}


inline Foam::tmp<Foam::volScalarField>
Foam::surfaceFilmModels::thermoSingleLayer::hs
(
    const volScalarField& T
) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "hs(" + T.name() + ")",
                filmRegion_.time().timeName(),
                filmRegion_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            cp_*(T - (dimensionedScalar("Tstd", dimTemperature, 298.15))),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}


inline Foam::tmp<Foam::volScalarField>
Foam::surfaceFilmModels::thermoSingleLayer::T
(
    const volScalarField& hs
) const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "T(" + hs.name() + ")",
                filmRegion_.time().timeName(),
                filmRegion_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            hs/cp_ + dimensionedScalar("Tstd", dimTemperature, 298.15),
            zeroGradientFvPatchScalarField::typeName
        )
    );
}


inline const Foam::volScalarField&
Foam::surfaceFilmModels::thermoSingleLayer::hsSp() const
{
    return hsSp_;
}


inline const Foam::volScalarField&
Foam::surfaceFilmModels::thermoSingleLayer::hsSpPrimary() const
{
    return hsSpPrimary_;
}


inline const Foam::volScalarField&
Foam::surfaceFilmModels::thermoSingleLayer::TPrimary() const
{
    return TPrimary_;
}


inline const Foam::PtrList<Foam::volScalarField>&
Foam::surfaceFilmModels::thermoSingleLayer::YPrimary() const
{
    return YPrimary_;
}


inline const Foam::surfaceFilmModels::heatTransferModel&
Foam::surfaceFilmModels::thermoSingleLayer::htcs() const
{
    return htcs_();
}


inline const Foam::surfaceFilmModels::heatTransferModel&
Foam::surfaceFilmModels::thermoSingleLayer::htcw() const
{
    return htcw_();
}


inline const Foam::surfaceFilmModels::phaseChangeModel&
Foam::surfaceFilmModels::thermoSingleLayer::phaseChange() const
{
    return phaseChange_();
}


// ************************************************************************* //
