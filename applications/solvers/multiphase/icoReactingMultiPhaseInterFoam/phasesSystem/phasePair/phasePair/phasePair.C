/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenCFD Ltd.
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

#include "phasePair.H"
#include "surfaceTensionModel.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::phasePair::phasePair
(
    const phaseModel& phase1,
    const phaseModel& phase2,
    const bool ordered
)
:
    phasePairKey(phase1.name(), phase2.name(), ordered),
    phase1_(phase1),
    phase2_(phase2),
    g_(phase1.mesh().lookupObject<uniformDimensionedVectorField>("g"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::phasePair::~phasePair()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::phaseModel& Foam::phasePair::dispersed() const
{
    FatalErrorIn("Foam::phasePair::dispersed() const")
        << "Requested dispersed phase from an unordered pair."
        << exit(FatalError);

    return phase1_;
}


const Foam::phaseModel& Foam::phasePair::continuous() const
{
    FatalErrorIn("Foam::phasePair::dispersed() const")
        << "Requested continuous phase from an unordered pair."
        << exit(FatalError);

    return phase1_;
}


Foam::word Foam::phasePair::name() const
{
    word name2(second());
    name2[0] = toupper(name2[0]);
    return first() + "And" + name2;
}


Foam::tmp<Foam::volScalarField> Foam::phasePair::rho() const
{
    return phase1()*phase1().rho() + phase2()*phase2().rho();
}


Foam::tmp<Foam::volScalarField> Foam::phasePair::Pr() const
{
    return
         continuous().nu()
        *continuous().Cp()
        *continuous().rho()
        /continuous().kappa();
}


Foam::tmp<Foam::volScalarField> Foam::phasePair::sigma() const
{
    return
        phase1().mesh().lookupObject<surfaceTensionModel>
        (
            IOobject::groupName
            (
                surfaceTensionModel::typeName,
                phasePair::name()
            )
        ).sigma();
}


// ************************************************************************* //
