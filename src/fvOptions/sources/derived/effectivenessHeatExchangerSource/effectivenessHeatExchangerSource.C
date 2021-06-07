/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2015 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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

#include "effectivenessHeatExchangerSource.H"
#include "fvMesh.H"
#include "fvMatrix.H"
#include "addToRunTimeSelectionTable.H"
#include "basicThermo.H"
#include "coupledPolyPatch.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(effectivenessHeatExchangerSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        effectivenessHeatExchangerSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fv::effectivenessHeatExchangerSource::initialise()
{
    const label zoneID = mesh_.faceZones().findZoneID(faceZoneName_);

    if (zoneID < 0)
    {
        FatalErrorInFunction
            << type() << " " << this->name() << ": "
            << "    Unknown face zone name: " << faceZoneName_
            << ". Valid face zones are: " << mesh_.faceZones().names()
            << exit(FatalError);
    }

    const faceZone& fZone = mesh_.faceZones()[zoneID];

    faceId_.setSize(fZone.size());
    facePatchId_.setSize(fZone.size());
    faceSign_.setSize(fZone.size());

    label count = 0;
    forAll(fZone, i)
    {
        label facei = fZone[i];
        label faceId = -1;
        label facePatchId = -1;
        if (mesh_.isInternalFace(facei))
        {
            faceId = facei;
            facePatchId = -1;
        }
        else
        {
            facePatchId = mesh_.boundaryMesh().whichPatch(facei);
            const polyPatch& pp = mesh_.boundaryMesh()[facePatchId];
            const auto* cpp = isA<coupledPolyPatch>(pp);

            if (cpp)
            {
                faceId = (cpp->owner() ? pp.whichFace(facei) : -1);
            }
            else if (!isA<emptyPolyPatch>(pp))
            {
                faceId = pp.whichFace(facei);
            }
            else
            {
                faceId = -1;
                facePatchId = -1;
            }
        }

        if (faceId >= 0)
        {
            if (fZone.flipMap()[i])
            {
                faceSign_[count] = -1;
            }
            else
            {
                faceSign_[count] = 1;
            }
            faceId_[count] = faceId;
            facePatchId_[count] = facePatchId;
            count++;
        }
    }
    faceId_.setSize(count);
    facePatchId_.setSize(count);
    faceSign_.setSize(count);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::effectivenessHeatExchangerSource::effectivenessHeatExchangerSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::cellSetOption(name, modelType, dict, mesh),
    secondaryMassFlowRate_(0),
    secondaryInletT_(0),
    primaryInletT_(0),
    userPrimaryInletT_(false),
    targetQdotActive_(false),
    targetQdot_(0),
    targetQdotCalcInterval_(5),
    targetQdotRelax_(0.5),
    eTable_(),
    UName_("U"),
    TName_("T"),
    phiName_("phi"),
    faceZoneName_("unknown-faceZone"),
    faceId_(),
    facePatchId_(),
    faceSign_()
{
    read(dict);

    // Set the field name to that of the energy
    // field from which the temperature is obtained

    const auto& thermo = mesh_.lookupObject<basicThermo>(basicThermo::dictName);

    fieldNames_.resize(1, thermo.he().name());

    fv::option::resetApplied();

    eTable_.reset(new interpolation2DTable<scalar>(coeffs_));

    initialise();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::effectivenessHeatExchangerSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label
)
{
    const auto& thermo = mesh_.lookupObject<basicThermo>(basicThermo::dictName);

    const surfaceScalarField Cpf(fvc::interpolate(thermo.Cp()));

    const auto& phi = mesh_.lookupObject<surfaceScalarField>(phiName_);

    const auto& T = mesh_.lookupObject<volScalarField>(TName_);
    const surfaceScalarField Tf(fvc::interpolate(T));

    scalar sumPhi = 0;
    scalar sumMagPhi = 0;
    scalar CpfMean = 0;
    scalar primaryInletTfMean = 0;
    forAll(faceId_, i)
    {
        label facei = faceId_[i];
        if (facePatchId_[i] != -1)
        {
            label patchi = facePatchId_[i];
            scalar phii = phi.boundaryField()[patchi][facei]*faceSign_[i];

            sumPhi += phii;

            scalar Cpfi = Cpf.boundaryField()[patchi][facei];
            scalar Tfi = Tf.boundaryField()[patchi][facei];
            scalar magPhii = mag(phii);

            sumMagPhi += magPhii;
            CpfMean += Cpfi*magPhii;
            primaryInletTfMean += Tfi*magPhii;
        }
        else
        {
            scalar phii = phi[facei]*faceSign_[i];
            scalar magPhii = mag(phii);

            sumPhi += phii;
            sumMagPhi += magPhii;
            CpfMean += Cpf[facei]*magPhii;
            primaryInletTfMean += Tf[facei]*magPhii;
        }
    }
    reduce(CpfMean, sumOp<scalar>());
    reduce(sumPhi, sumOp<scalar>());
    reduce(sumMagPhi, sumOp<scalar>());
    CpfMean /= sumMagPhi + ROOTVSMALL;

    scalar primaryInletT = primaryInletT_;
    if (!userPrimaryInletT_)
    {
        reduce(primaryInletTfMean, sumOp<scalar>());
        primaryInletT = primaryInletTfMean/(sumMagPhi + ROOTVSMALL);
    }

    const scalar alpha =
        eTable_()(mag(sumPhi), secondaryMassFlowRate_)
       *CpfMean*mag(sumPhi);

    const scalar Qt = alpha*(secondaryInletT_ - primaryInletT);

    if
    (
        targetQdotActive_
     && (mesh_.time().timeIndex() % targetQdotCalcInterval_ == 0)
    )
    {
        scalar dT = (targetQdot_ - Qt)/(alpha + ROOTVSMALL);
        secondaryInletT_ += targetQdotRelax_*dT;
    }

    const scalarField TCells(T, cells_);
    scalar Tref = 0;
    scalarField deltaTCells(cells_.size(), Zero);
    if (Qt > 0)
    {
        Tref = gMax(TCells);
        forAll(deltaTCells, i)
        {
            deltaTCells[i] = max(Tref - TCells[i], 0.0);
        }
    }
    else
    {
        Tref = gMin(TCells);
        forAll(deltaTCells, i)
        {
            deltaTCells[i] = max(TCells[i] - Tref, 0.0);
        }
    }


    const auto& U = mesh_.lookupObject<volVectorField>(UName_);
    const scalarField& V = mesh_.V();
    scalar sumWeight = 0;
    forAll(cells_, i)
    {
        label celli = cells_[i];
        sumWeight += V[celli]*mag(U[celli])*deltaTCells[i];
    }
    reduce(sumWeight, sumOp<scalar>());

    if (this->V() > VSMALL && mag(Qt) > VSMALL)
    {
        scalarField& heSource = eqn.source();

        forAll(cells_, i)
        {
            label celli = cells_[i];
            heSource[celli] -=
                Qt*V[celli]*mag(U[celli])*deltaTCells[i]
               /(sumWeight + ROOTVSMALL);
        }
    }

    Info<< type() << ": " << name() << nl << incrIndent
        << indent << "Net mass flux [Kg/s]      : " << sumPhi << nl
        << indent << "Total heat exchange [W] : " << Qt << nl
        << indent << "Secondary inlet T [K]     : " << secondaryInletT_ << nl
        << indent << "Tref [K]                  : " << Tref << nl
        << indent << "Effectiveness             : "
        << eTable_()(mag(sumPhi), secondaryMassFlowRate_) << decrIndent
        << nl << endl;
}


bool Foam::fv::effectivenessHeatExchangerSource::read(const dictionary& dict)
{
    if (fv::cellSetOption::read(dict))
    {
        UName_ = coeffs_.getOrDefault<word>("U", "U");
        TName_ = coeffs_.getOrDefault<word>("T", "T");
        phiName_ = coeffs_.getOrDefault<word>("phi", "phi");
        coeffs_.readEntry("faceZone", faceZoneName_);

        coeffs_.readEntry("secondaryMassFlowRate", secondaryMassFlowRate_);
        coeffs_.readEntry("secondaryInletT", secondaryInletT_);

        if (coeffs_.readIfPresent("primaryInletT", primaryInletT_))
        {
            userPrimaryInletT_ = true;
            Info<< type() << " " << this->name() << ": " << indent << nl
                << "employing user-specified primary flow inlet temperature: "
                << primaryInletT_ << endl;
        }
        else
        {
            Info<< type() << " " << this->name() << ": " << indent << nl
                << "employing flux-weighted primary flow inlet temperature"
                << endl;
        }

        if (coeffs_.readIfPresent("targetQdot", targetQdot_))
        {
            targetQdotActive_ = true;
            Info<< indent << "employing target heat rejection of "
                << targetQdot_ << nl;

            coeffs_.readIfPresent
            (
                "targetQdotCalcInterval",
                targetQdotCalcInterval_
            );

            Info<< indent << "updating secondary inlet temperature every "
                << targetQdotCalcInterval_ << " iterations" << nl;

            coeffs_.readIfPresent("targetQdotRelax", targetQdotRelax_);

            Info<< indent << "temperature relaxation:  "
                << targetQdotRelax_ << endl;
        }

        return true;
    }

    return false;
}


// ************************************************************************* //
