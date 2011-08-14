/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::Time& Foam::spray::runTime() const
{
    return runTime_;
}


inline const Foam::fvMesh& Foam::spray::mesh() const
{
    return mesh_;
}


inline const Foam::volVectorField& Foam::spray::U() const
{
    return U_;
}


inline const Foam::volScalarField& Foam::spray::rho() const
{
    return rho_;
}


inline const Foam::volScalarField& Foam::spray::p() const
{
    return p_;
}


inline const Foam::volScalarField& Foam::spray::T() const
{
    return T_;
}


inline Foam::PtrList<Foam::injector>& Foam::spray::injectors()
{
    return injectors_;
}


inline const Foam::PtrList<Foam::injector>& Foam::spray::injectors() const
{
    return injectors_;
}


inline const Foam::atomizationModel& Foam::spray::atomization() const
{
    return atomization_();
}


inline const Foam::breakupModel& Foam::spray::breakup() const
{
    return breakupModel_();
}


inline const Foam::collisionModel& Foam::spray::collisions() const
{
    return collisionModel_();
}


inline const Foam::dispersionModel& Foam::spray::dispersion() const
{
    return dispersionModel_();
}


inline const Foam::dragModel& Foam::spray::drag() const
{
    return drag_();
}


inline const Foam::evaporationModel& Foam::spray::evaporation() const
{
    return evaporation_();
}


inline const Foam::heatTransferModel& Foam::spray::heatTransfer() const
{
    return heatTransfer_();
}


inline const Foam::injectorModel& Foam::spray::injection() const
{
    return injectorModel_();
}


inline const Foam::wallModel& Foam::spray::wall() const
{
    return wall_();
}


inline Foam::tmp<Foam::volVectorField> Foam::spray::momentumSource() const
{
    tmp<volVectorField> tsource
    (
        new volVectorField
        (
            IOobject
            (
                "sms",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector
            (
                "zero",
                dimensionSet(1, -2, -2, 0, 0),
                vector::zero
            )
        )
    );

    tsource().internalField() = sms_/runTime_.deltaTValue()/mesh_.V();

    return tsource;
}


inline Foam::tmp<Foam::volScalarField>
Foam::spray::evaporationSource(const label si) const
{
    tmp<volScalarField> tsource
    (
        new volScalarField
        (
            IOobject
            (
                "srhos",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimensionSet(1, -3, -1, 0, 0), 0.0)
        )
    );

    if (isLiquidFuel_[si])
    {
        label fi = gasToLiquidIndex_[si];
        tsource().internalField() = srhos_[fi]/runTime_.deltaTValue()/mesh_.V();
    }
    else
    {
        scalarField s(mesh_.nCells(), 0.0);
        tsource().internalField() = s;
    }

    return tsource;
}


inline Foam::tmp<Foam::volScalarField> Foam::spray::heatTransferSource() const
{
    tmp<volScalarField> tsource
    (
        new volScalarField
        (
            IOobject
            (
                "shs",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimensionSet(1, -1, -3, 0, 0), 0.0)
        )
    );

    tsource().internalField() = shs_/runTime_.deltaTValue()/mesh_.V();

    return tsource;
}


inline Foam::cachedRandom& Foam::spray::rndGen()
{
    return rndGen_;
}


inline Foam::label Foam::spray::subCycles() const
{
    return subCycles_;
}


inline const Foam::vector& Foam::spray::g() const
{
    return g_;
}


inline const Foam::liquidMixtureProperties& Foam::spray::fuels() const
{
    return fuels_();
}


inline const Foam::PtrList<Foam::gasThermoPhysics>&
Foam::spray::gasProperties() const
{
    return gasProperties_;
}


inline const Foam::basicMultiComponentMixture& Foam::spray::composition() const
{
    return composition_;
}


inline const Foam::List<Foam::label>& Foam::spray::liquidToGasIndex() const
{
    return liquidToGasIndex_;
}


inline const Foam::List<Foam::label>& Foam::spray::gasToLiquidIndex() const
{
    return gasToLiquidIndex_;
}


inline const Foam::List<bool>& Foam::spray::isLiquidFuel() const
{
    return isLiquidFuel_;
}


inline bool Foam::spray::twoD() const
{
    return twoD_;
}


inline const Foam::vector& Foam::spray::axisOfSymmetry() const
{
    return axisOfSymmetry_;
}


inline const Foam::vector& Foam::spray::axisOfWedge() const
{
    return axisOfWedge_;
}


inline const Foam::vector& Foam::spray::axisOfWedgeNormal() const
{
    return axisOfWedgeNormal_;
}


inline Foam::scalar Foam::spray::angleOfWedge() const
{
    return angleOfWedge_;
}


inline const Foam::interpolation<Foam::vector>&
Foam::spray::UInterpolator() const
{
    return UInterpolator_;
}


inline const Foam::interpolation<Foam::scalar>&
Foam::spray::rhoInterpolator() const
{
    return rhoInterpolator_;
}


inline const Foam::interpolation<Foam::scalar>&
Foam::spray::pInterpolator() const
{
    return pInterpolator_;
}


inline const Foam::interpolation<Foam::scalar>&
Foam::spray::TInterpolator() const
{
    return TInterpolator_;
}


inline Foam::vectorField& Foam::spray::sms()
{
    return sms_;
}


inline const Foam::vectorField& Foam::spray::sms() const
{
    return sms_;
}


inline Foam::scalarField& Foam::spray::shs()
{
    return shs_;
}


inline const Foam::scalarField& Foam::spray::shs() const
{
    return shs_;
}


inline Foam::PtrList<Foam::scalarField>& Foam::spray::srhos()
{
    return srhos_;
}


inline const Foam::PtrList<Foam::scalarField>& Foam::spray::srhos() const
{
    return srhos_;
}


inline Foam::scalar Foam::spray::ambientPressure() const
{
    return ambientPressure_;
}


inline Foam::scalar Foam::spray::ambientTemperature() const
{
    return ambientTemperature_;
}


// ************************************************************************* //
