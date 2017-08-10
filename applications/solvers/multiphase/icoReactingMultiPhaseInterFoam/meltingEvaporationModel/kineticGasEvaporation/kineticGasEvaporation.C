/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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

#include "kineticGasEvaporation.H"
#include "constants.H"
#include "fvcGrad.H"
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
    C_("C",  dimless, dict.lookup("C")),
    Tactivate_
    (
        "Tactivate",
        dimTemperature,
        dict.lookup("Tactivate")
    ),
    Mv_
    (
        "Mv",
        dimMass/dimMoles,
        dict.lookupOrDefault<scalar>("Mv", -1)
    ),
    saturationPressureModel_
    (
        saturationPressureModel::New
        (
            dict.subDict("saturationPressure")
        )
    ),
    phi0_
    (
        IOobject
        (
            "phi0",
            this->mesh_.time().timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("zero", dimMass/dimTime/dimVolume, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    sPhi_
    (
        IOobject
        (
            "sPhi",
            this->mesh_.time().timeName(),
            this->mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("zero", dimMass/dimTime/dimVolume, 0),
        zeroGradientFvPatchScalarField::typeName
    ),
    D_
    (
        "D",
        dimArea/dimTime,
        dict.lookupOrDefault<scalar>("D", 1)
    )
{
    word speciesName = this->species()[0];

    // Get the continuous thermo
    const typename OtherThermo::thermoType& localThermo =
        this->getLocalThermo
        (
            speciesName,
            this->otherThermo_
        );

    Mv_.value() = localThermo.W();

    if (Mv_.value() == -1)
    {
        FatalErrorIn
        (
           "meltingEvaporationModels::"
           "kineticGasEvaporation<Thermo, OtherThermo>::"
           "kineticGasEvaporation"
           "("
            "const dictionary& dict,"
            "const phasePair& pair"
           ")"
        )
            << " Please provide the molar weight (Mv) of vapour [g/mol] "
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::meltingEvaporationModels::kineticGasEvaporation<Thermo, OtherThermo>
::Kexp(label variable, const volScalarField& field) const
{

    if (this->modelVariable_ == variable)
    {
        volScalarField limitedContinous
        (
            min(max(this->pair().continuous(), scalar(0)), scalar(1))
        );

        volScalarField limitedDispersed
        (
            min(max(this->pair().dispersed(), scalar(0)), scalar(1))
        );

        const fvMesh& mesh = this->mesh_;


        const volScalarField& T = mesh.lookupObject<volScalarField>("T");
        //const volScalarField& p = mesh.lookupObject<volScalarField>("p");

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
/*
        volScalarField pSat
        (
            "pSat",
            saturationPressureModel_->pSat(T)
        );
*/
/*
        volScalarField rhoMean
        (
            this->pair().continuous().rho()*this->pair().dispersed().rho()
        /  (this->pair().dispersed().rho() - this->pair().continuous().rho())
        );
*/
        word species(this->species()[0]);

        tmp<volScalarField> L = mag(this->L(species, field));

        DebugVar(Mv_);
        //dimensionedScalar psat("psat", dimPressure, 1e5);
        //scalarField maxAreaDen(sqrt(3.0)*pow(mesh.V(), -1.0/3.0));
        //scalarField avergageDelta(pow(mesh.V(), 1.0/3.0));
/*
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
*/
        volScalarField areaDensity
        (
            "areaDensity",
            mag(fvc::grad(limitedDispersed))
        );

/*
        volScalarField tempInterFace
        (
            "tempInterFace",
            pos(0.55 - limitedDispersed) * pos(limitedDispersed - 0.45)
          * pos(limitedContinous - 0.45) * pos(0.55 - limitedContinous)
          * T
        );
*/
        dimensionedScalar dMgc = max(areaDensity);


        volScalarField Tave
        (
            "Tave",  T*pos(areaDensity - 0.1*dMgc)
        );

        const volScalarField massFluxEvap
        (
            "massFluxEvap",
            2*C_/(2 - C_)
          * HerztKnudsConst
          * L()
          * this->pair().continuous().rho()
          * max
            (
                (Tave - Tactivate_),
                dimensionedScalar("T0", dimTemperature, 0.0)
            )
        );

        //scalar alphaC(0.9);
        //scalar alphaCMax(0.9);


        //.weightedAverage(mesh.V());
/*
        volScalarField nearInterFace
        (
            "nearInterFace",
            pos(areaDensity - 0.1*dMgc)
          * pos(alphaC - limitedDispersed)
          * pos(limitedDispersed - (1 - alphaC))
          * pos(limitedContinous - (1 - alphaC))
          * pos(alphaC - limitedContinous)
        );

        volScalarField farDispersedInterFace
        (
            pos(limitedDispersed - alphaC)
        );

        volScalarField farContinousInterFace
        (
            pos(limitedContinous - alphaC)
        );

        volScalarField farInterFace
        (
            "farInterFace", farDispersedInterFace + farContinousInterFace
        );
*/
        dimensionedScalar mIntDotPhi0
        (
            "mIntDotPhi0",
            gSum((pos(areaDensity - 0.1*dMgc)*areaDensity*mesh.V())())
           /gSum
           (
                (pos(areaDensity - 0.1*dMgc)*limitedDispersed*areaDensity*mesh.V())()
           )
        );

        DebugVar(mIntDotPhi0);

        // Local density rate (kg/m3/s)
        phi0_ = massFluxEvap*areaDensity; //mIntDotPhi0*limitedDispersed


        dimensionedScalar mIntDot("mIntDot", gSum((phi0_*mesh.V())()));

        DebugVar(mIntDot);

        volScalarField sPhi("sPhi", phi0_);

        phi0_ = 0.0*phi0_;

        fvScalarMatrix sPhiEq
        (
            fvm::ddt(sPhi)
          - fvm::laplacian(D_, sPhi)
        );

        sPhiEq.relax();
        sPhiEq.solve();

        // Liquid normalization
        const dimensionedScalar Nl
        (
            gSum((sPhi*mesh.V())())
           /(
               gSum
               (
                   (sPhi*mesh.V()*limitedDispersed)()
               )
             + dimensionedScalar("SMALL", dimless, VSMALL)
            )
        );

        // Vapour normalization
        const dimensionedScalar Nv
        (
            gSum((sPhi*mesh.V())())
           /(
               gSum
               (
                   (sPhi*mesh.V()*limitedContinous)()
               )
             + dimensionedScalar("SMALL", dimless, VSMALL)
            )
        );

        // Spread density source
        sPhi_ =
            0.5*
            (
                limitedContinous*Nv*sPhi
              + limitedDispersed*Nl*sPhi
            );

        dimensionedScalar msPhiIntDot("msPhiIntDot", gSum((sPhi_*mesh.V())()));

        DebugVar(msPhiIntDot);

        if (this->pair().dispersed().mesh().time().outputTime())
        {
            massFluxEvap.write();
            areaDensity.write();
            Tave.write();
        }

        return
        (
            0.5*
            (
                limitedContinous*Nv*sPhi
              + limitedDispersed*Nl*sPhi
            )
        );
    }
    else
    {
        return tmp<volScalarField> ();
    }
}


template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::meltingEvaporationModels::kineticGasEvaporation<Thermo, OtherThermo>
::Kimp(label variable, const volScalarField& field) const
{
    if (this->modelVariable_ == variable)
    {
         return tmp<volScalarField> ();

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


template<class Thermo, class OtherThermo>
Foam::tmp<Foam::volScalarField>
Foam::meltingEvaporationModels::kineticGasEvaporation<Thermo, OtherThermo>
::KexpEnergy
(
    label modelVariable,
    const volScalarField& field
)
const
{
    return dimensionedScalar("one", dimless, 1.0)*phi0_;
}


template<class Thermo, class OtherThermo>
Foam::label
Foam::meltingEvaporationModels::kineticGasEvaporation<Thermo, OtherThermo>
::dSdVariable()
{
    return label(1);
}

// ************************************************************************* //
