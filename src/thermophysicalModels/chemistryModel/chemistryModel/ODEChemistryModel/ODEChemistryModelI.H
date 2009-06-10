/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "zeroGradientFvPatchFields.H"

template<class ThermoType>
inline Foam::PtrList<Foam::scalarField>&
Foam::ODEChemistryModel<ThermoType>::RR()
{
    return RR_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoType>
inline const Foam::PtrList<Foam::Reaction<ThermoType> >&
Foam::ODEChemistryModel<ThermoType>::reactions() const
{
    return reactions_;
}


template<class ThermoType>
inline const Foam::PtrList<ThermoType>&
Foam::ODEChemistryModel<ThermoType>::specieThermo() const
{
    return specieThermo_;
}


template<class ThermoType>
inline Foam::label Foam::ODEChemistryModel<ThermoType>::nSpecie() const
{
    return nSpecie_;
}


template<class ThermoType>
inline Foam::label Foam::ODEChemistryModel<ThermoType>::nReaction() const
{
    return nReaction_;
}


template<class ThermoType>
inline const Foam::chemistrySolver<ThermoType>&
Foam::ODEChemistryModel<ThermoType>::solver() const
{
    return solver_;
}


template<class ThermoType>
inline Foam::tmp<Foam::volScalarField>
Foam::ODEChemistryModel<ThermoType>::RR
(
    const label i
) const
{
    tmp<volScalarField> tRR
    (
        new volScalarField
        (
            IOobject
            (
                "RR(" + this->Y_[i].name() + ')',
                this->thermo().rho()().mesh().time().timeName(),
                this->thermo().rho()().mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            this->thermo().rho()().mesh(),
            dimensionedScalar("zero", dimensionSet(1, -3, -1, 0, 0), 0.0),
            zeroGradientFvPatchScalarField::typeName
        )
    );

    if (chemistry_)
    {
        tRR().internalField() = RR_[i];
        tRR().correctBoundaryConditions();
    }
    return tRR;
}


// ************************************************************************* //
