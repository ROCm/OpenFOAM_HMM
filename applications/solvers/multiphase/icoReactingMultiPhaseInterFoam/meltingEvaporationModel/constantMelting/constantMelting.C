/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd
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

#include "constantMelting.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::meltingEvaporationModels::constantMelting<Thermo, OtherThermo>::
constantMelting
(
    const dictionary& dict,
    const phasePair& pair
)
:
    InterfaceCompositionModel<Thermo, OtherThermo>(dict, pair),
    C_("C",  inv(dimTime)*inv(dimTemperature), dict.lookup("C")),
    Tactivate_("Tactivate", dimTemperature, dict.lookup("Tactivate"))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::meltingEvaporationModels::constantMelting<Thermo, OtherThermo>
::Kexp(label variable, const volScalarField& field) const
{

    if (this->modelVariable_ == variable)
    {
        volScalarField limitedDispersed
        (
            min(max(this->pair().dispersed(), scalar(0)), scalar(1))
        );

        return
        (
            C_
          * limitedDispersed
          * this->pair().dispersed().rho()
          * (field.oldTime() - Tactivate_)
          * pos(field.oldTime() - Tactivate_)
        );
    }
    else
    {
        return tmp<volScalarField> ();
    }
}

template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::meltingEvaporationModels::constantMelting<Thermo, OtherThermo>
::Kimp(label variable, const volScalarField& field) const
{
    if (this->modelVariable_ == variable)
    {
        volScalarField limitedDispersed
        (
            min(max(this->pair().dispersed(), scalar(0)), scalar(1))
        );

        return
        (
            C_
          * limitedDispersed
          * this->pair().dispersed().rho()
          * pos(field.oldTime() - Tactivate_)
        );
    }
    else
    {
        return tmp<volScalarField> ();
    }

}


template<class Thermo, class OtherThermo>
const Foam::dimensionedScalar&
Foam::meltingEvaporationModels::constantMelting<Thermo, OtherThermo>
::Tactivate() const
{
    return Tactivate_;
}


template<class Thermo, class OtherThermo>
Foam::label
Foam::meltingEvaporationModels::constantMelting<Thermo, OtherThermo>
::dSdVariable()
{
    return label(1);
}

// ************************************************************************* //
