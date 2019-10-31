/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenCFD Ltd.
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

#include "kineticGasEvaporation.H"
#include "constants.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcDiv.H"
#include "surfaceInterpolate.H"
#include "fvcReconstruct.H"
#include "fvm.H"
#include "zeroGradientFvPatchFields.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::meltingEvaporationModels::kineticGasEvaporation<Thermo, OtherThermo>
::kineticGasEvaporation
(
    const dictionary& dict,
    const phasePair& pair
)
:
    InterfaceCompositionModel<Thermo, OtherThermo>(dict, pair),
    C_("C", dimless, dict),
    Tactivate_("Tactivate", dimTemperature, dict),
    Mv_("Mv", dimMass/dimMoles, -1, dict),
    alphaMax_(dict.lookupOrDefault<scalar>("alphaMax", 1.0)),
    alphaMin_(dict.lookupOrDefault<scalar>("alphaMin", 0.5)),
    alphaRestMax_(dict.lookupOrDefault<scalar>("alphaRestMax", 0.01))
{
    if (this->transferSpecie() != "none")
    {
        word fullSpeciesName = this->transferSpecie();
        auto tempOpen = fullSpeciesName.find('.');
        const word speciesName(fullSpeciesName.substr(0, tempOpen));

        // Get the "to" thermo
        const typename OtherThermo::thermoType& toThermo =
            this->getLocalThermo
            (
                speciesName,
                this->toThermo_
            );

         // Convert from g/mol to Kg/mol
        Mv_.value() = toThermo.W()*1e-3;
    }


    if (Mv_.value() == -1)
    {
        FatalErrorInFunction
            << " Please provide the molar weight (Mv) of vapour [g/mol] "
            << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::meltingEvaporationModels::kineticGasEvaporation<Thermo, OtherThermo>
::Kexp(label variable, const volScalarField& field)
{
    if (this->modelVariable_ == variable)
    {
        const volScalarField& to = this->pair().to();

        const volScalarField& from = this->pair().from();

        const fvMesh& mesh = this->mesh_;

        const volScalarField& T =
            mesh.lookupObject<volScalarField>("T").oldTime();

        const dimensionedScalar HerztKnudsConst
        (
            sqrt
            (
                Mv_
               /2.0
               /constant::physicoChemical::R
               /mathematical::pi
               /pow3(Tactivate_)
            )
        );

        word fullSpeciesName = this->transferSpecie();
        auto tempOpen = fullSpeciesName.find('.');
        const word speciesName(fullSpeciesName.substr(0, tempOpen));

        tmp<volScalarField> L = this->L(speciesName, field);

        const volVectorField gradFrom(fvc::grad(from));
        const volVectorField gradTo(fvc::grad(to));

        const volScalarField areaDensity("areaDensity", mag(gradFrom));

        const volScalarField gradAlphaf(gradFrom & gradTo);

        volScalarField Tmask("Tmask", from*0.0);

        forAll(Tmask, celli)
        {
            if (gradAlphaf[celli] < 0)
            {
                if (from[celli] > alphaMin_ && from[celli] < alphaMax_)
                {
                    {
                        scalar alphaRes = 1.0 - from[celli] - to[celli];
                        if (alphaRes < alphaRestMax_)
                        {
                            Tmask[celli] = 1.0;
                        }
                    }
                }
            }
        }

        tmp<volScalarField> tRhom
        (
            new volScalarField
            (
                IOobject
                (
                    "trhom",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar(dimDensity, Zero)
            )
        );
        volScalarField& rhom = tRhom.ref();

        tmp<volScalarField> tTdelta
        (
            new volScalarField
            (
                IOobject
                (
                    "trhom",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar(dimTemperature, Zero)
            )
        );
        volScalarField& tDelta = tTdelta.ref();

        if (sign(C_.value()) > 0)
        {
            rhom =
                this->pair().to().rho()*this->pair().from().rho()
              / (this->pair().from().rho() - this->pair().to().rho());

            tDelta = max
            (
                (T*Tmask - Tactivate_),
                dimensionedScalar("T0", dimTemperature, Zero)
            );
        }
        else
        {
            rhom =
                this->pair().to().rho()*this->pair().from().rho()
              / (this->pair().to().rho() - this->pair().from().rho());

            tDelta = max
            (
                Tmask*(Tactivate_ - T),
                dimensionedScalar("T0", dimTemperature, Zero)
            );
        }

        volScalarField massFluxEvap
        (
            "massFluxEvap",
            2*mag(C_)/(2 - mag(C_))
          * HerztKnudsConst
          * L()
          * rhom
          * tDelta
        );

        // 'from' phase normalization
        // WIP: Normalization could be convinient for cases where the area were
        // the source term is calculated is uniform
        const dimensionedScalar Nl
        (
            gSum((areaDensity*mesh.V())())
           /(
               gSum
               (
                   ((areaDensity*from)*mesh.V())()
               )
             + dimensionedScalar("SMALL", dimless, VSMALL)
            )
        );


        if (mesh.time().outputTime() && debug)
        {
            areaDensity.write();
            Tmask.write();
            volScalarField mKGasDot
            (
                "mKGasDot",
                massFluxEvap*areaDensity*Nl*from
            );
            mKGasDot.write();
        }

        return massFluxEvap*areaDensity*Nl*from;
    }
    else
    {
        return tmp<volScalarField> ();
    }
}


template<class Thermo, class OtherThermo>
const Foam::dimensionedScalar&
Foam::meltingEvaporationModels::kineticGasEvaporation<Thermo, OtherThermo>
::Tactivate() const
{
    return Tactivate_;
}


// ************************************************************************* //
