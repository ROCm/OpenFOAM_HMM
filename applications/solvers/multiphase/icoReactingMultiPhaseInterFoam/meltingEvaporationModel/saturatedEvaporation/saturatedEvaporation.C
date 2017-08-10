/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C)  2015 OpenCFD Ltd.
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

#include "saturatedEvaporation.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::meltingEvaporationModels::saturatedEvaporation<Thermo, OtherThermo>::
saturatedEvaporation
(
    const dictionary& dict,
    const phasePair& pair
)
:
    InterfaceCompositionModel<Thermo, OtherThermo>(dict, pair),
    C_("C",  inv(dimTime), dict.lookup("C")),
    Tactivate_("Tactivate", dimTemperature, dict.lookup("Tactivate")),
    saturationPressureModel_
    (
        saturationPressureModel::New
        (
            dict.subDict("saturationPressure")
        )
    )
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField> Foam::meltingEvaporationModels::
saturatedEvaporation<Thermo, OtherThermo>::Yf
(
    const word& speciesName,
    const volScalarField& Tf
) const
{
    const fvMesh& mesh = this->pair().dispersed().mesh();

    tmp<volScalarField> twRatio
    (
        new volScalarField
        (
            IOobject
            (
                "wRatio",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("wRatio", dimless, 1)
        )
    );

    volScalarField& wRatio = twRatio();

    // Get the continuous thermo
    const typename OtherThermo::thermoType& continuousThermo =
        this->getLocalThermo
        (
            speciesName,
            this->otherThermo_
        );

//     // Get the dispersed thermo
//     const typename Thermo::thermoType& dispersedThermo =
//         this->getLocalThermo
//         (
//             speciesName,
//             this->thermo_
//         );

    // If dispersed phase (liquid) is multicomponent
    {
        wRatio =
            this->MwMixture(this->thermo_)
          * this->getSpecieMassFraction(speciesName, this->thermo_)
          / this->MwMixture(this->thermo_);
    }

    // If continuous phase (vapour) is multicomponent
    {
        wRatio *=
            continuousThermo.W()
          / this->MwMixture(this->otherThermo_);
    }

    const volScalarField& p = this->pair().dispersed().thermo().p();

    volScalarField Yf
    (
        "Yf",
        wRatio*saturationPressureModel_->pSat(Tf)/p
    );

    volScalarField pSat
    (
        "pSat",
         saturationPressureModel_->pSat(Tf)
    );

    if (this->pair().dispersed().mesh().time().outputTime())
    {
        Yf.write();
        pSat.write();
    }

    return wRatio*saturationPressureModel_->pSat(Tf)/p;
}


template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::meltingEvaporationModels::saturatedEvaporation<Thermo, OtherThermo>
::Kexp(label variable, const volScalarField& field) const
{
    if (this->modelVariable_ == variable)
    {
        volScalarField limitedDispersed
        (
            min(max(this->pair().dispersed(), scalar(0)), scalar(1))
        );

        volScalarField limitedContinous
        (
            min(max(this->pair().continuous(), scalar(0)), scalar(1))
        );

        volScalarField nearInterFace
        (
            "nearInterFace",
            pos(limitedDispersed - 0.1)*pos(0.9 - limitedDispersed)
           *pos(limitedContinous - 0.1)*pos(0.9 - limitedContinous)
        );

        return
        (
            C_
           *nearInterFace
           *limitedDispersed
           *pos(field.oldTime() - Tactivate_)
           *this->pair().dispersed().rho()
        );
    }
    else
    {
         return tmp<volScalarField> ();
    }
}


template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::meltingEvaporationModels::saturatedEvaporation<Thermo, OtherThermo>
::Kimp(label variable, const volScalarField& field) const
{
    if (this->modelVariable_ == variable)
    {
        /*
        volScalarField limitedDispersed
        (
            min(max(this->pair().dispersed(), scalar(0)), scalar(1))
        );

        volScalarField limitedContinous
        (
            min(max(this->pair().continuous(), scalar(0)), scalar(1))
        );

        volScalarField nearInterFace
        (
            "nearInterFace",
            pos(limitedDispersed - 0.1)*pos(0.9 - limitedDispersed)
           *pos(limitedContinous - 0.1)*pos(0.9 - limitedContinous)
        );

        scalarField maxAreaDen(sqrt(3.0)*pow(mesh.V(), -1.0/3.0));

        volScalarField areaDensity
        (
            IOobject
            (
                "areaDensity",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("zero", dimArea/dimVolume, 0.0),
            zeroGradientFvPatchScalarField::typeName
        );

        forAll (limitedDispersed, i)
        {
            if (limitedDispersed[i] < 0.5)
            {
                areaDensity[i] = 2.0*limitedDispersed[i]*maxAreaDen[i];
            }
            else
            {
                areaDensity[i] = 2.0*(1.0 - limitedDispersed[i])*maxAreaDen[i];
            }
        }
        areaDensity.correctBoundaryConditions();

        return
        (
            C_
           *nearInterFace
           *limitedDispersed
           *areaDensity
           *pos(field.oldTime() - Tactivate_)
        );
        */
         return tmp<volScalarField> ();
    }
    else
    {
        return tmp<volScalarField> ();
    }
}


template<class Thermo, class OtherThermo>
const Foam::dimensionedScalar&
Foam::meltingEvaporationModels::saturatedEvaporation<Thermo, OtherThermo>
::Tactivate() const
{
    return Tactivate_;
}

// ************************************************************************* //
