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

#include "interRegionHeatTransferModel.H"
#include "basicThermo.H"
#include "fvm.H"
#include "IObasicSourceList.H"
#include "zeroGradientFvPatchFields.H"
#include "fvcVolumeIntegrate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(interRegionHeatTransferModel, 0);
};


// * * * * * * * * * * * *  Private member functions * * * * * * * * * * * //

void Foam::interRegionHeatTransferModel::check()
{
    const fvMesh& secondaryMesh =
            mesh_.time().lookupObject<fvMesh>(mapRegionName_);

    const basicSourceList& IObsl =
        secondaryMesh.lookupObject<basicSourceList>
        (
            "sourcesProperties"
        );

    const PtrList<basicSource>& bsl = IObsl;

    bool secSourceFound(false);

    forAll(bsl, i)
    {
        if (bsl[i].name() == secondarySourceName_)
        {
            secIrht_ = &const_cast<interRegionHeatTransferModel&>
            (
                refCast<const interRegionHeatTransferModel>(bsl[i])
            );
            secSourceFound = true;
            break;
        }
    }

    if(!secSourceFound)
    {
        FatalErrorIn
        (
            "constantHeatTransfer::interRegionHeatTransferModel"
            "("
            "   const word& name,"
            "   const word& modelType,"
            "   const dictionary& dict,"
            "   const fvMesh& mesh"
            ")"
        )   << "Secondary source name not found" << secondarySourceName_
            << " in region " << secondaryMesh.name()
            << nl
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interRegionHeatTransferModel::interRegionHeatTransferModel
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    basicSource(name, modelType, dict, mesh),
    secIrht_(),
    firstIter_(true),
    htc_
    (
        IOobject
        (
            "htc",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar
        (
            "htc",
            dimEnergy/dimTime/dimTemperature/dimVolume,
            0.0
        ),
        zeroGradientFvPatchScalarField::typeName
    )
{
    coeffs_.lookup("fieldNames") >> fieldNames_;
    applied_.setSize(fieldNames_.size(), false);
}


Foam::interRegionHeatTransferModel::~interRegionHeatTransferModel()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::interRegionHeatTransferModel::addSup
(
    fvMatrix<scalar>& eEqn,
    const label fieldI
)
{
    if (firstIter_)
    {
        check();
        firstIter_ = false;
    }

    const volScalarField& h = eEqn.psi();

    tmp<volScalarField> tTmapped
    (
        new volScalarField
        (
            IOobject
            (
                "Tmapped" + mesh_.name(),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("T", dimTemperature, 0.0)
        )
    );

    volScalarField& Tmapped = tTmapped();

    const fvMesh& secondaryMesh =
        mesh_.time().lookupObject<fvMesh>(mapRegionName_);

    basicThermo&  secondaryThermo =
        const_cast<basicThermo&>
        (
            secondaryMesh.lookupObject<basicThermo>("thermophysicalProperties")
        );

    secondaryToPrimaryInterpPtr_->interpolateInternalField
    (
        Tmapped,
        secondaryThermo.T(),
        meshToMesh::MAP,
        eqOp<scalar>()
    );

    if (!master_)
    {
        secondaryToPrimaryInterpPtr_->interpolateInternalField
        (
            htc_,
            secIrht_->calculateHtc(),
            meshToMesh::CELL_VOLUME_WEIGHT,
            eqOp<scalar>()
        );
    }

    if (debug)
    {
        Info<< " Volumetric integral of htc : "
            << fvc::domainIntegrate(htc_).value()
            << endl;
    }

    //SAF: temporarily output
    if (mesh_.time().outputTime())
    {
        Tmapped.write();
        htc_.write();
    }

    if (h.dimensions() == dimEnergy/dimMass)
    {
        const basicThermo& primaryThermo =
            mesh_.lookupObject<basicThermo>("thermophysicalProperties");

        eEqn += htc_*Tmapped - fvm::Sp(htc_/primaryThermo.Cp(), h);

        if (debug)
        {
            Info<< " Energy exchange from region " << secondaryMesh.name()
                << " To " << mesh_.name() << " : "
                <<  fvc::domainIntegrate
                    (
                        htc_*(h/primaryThermo.Cp() - Tmapped)
                    ).value()
                << endl;
        }
    }
    else if(h.dimensions() == dimTemperature)
    {
        eEqn += htc_*Tmapped - fvm::Sp(htc_, h);

        if (debug)
        {
            Info<< " Enegy exchange from region " << secondaryMesh.name()
                << " To " << mesh_.name() << " : "
                <<  fvc::domainIntegrate(htc_*(h - Tmapped)).value()
                << endl;
        }
    }
}


void Foam::interRegionHeatTransferModel::writeData(Ostream& os) const
{
    os.writeKeyword("name") << this->name() << token::END_STATEMENT << nl;
    os.writeKeyword("mapRegionName") << mapRegionName_
        << token::END_STATEMENT << nl;
    os.writeKeyword("secondarySourceName") << secondarySourceName_
        << token::END_STATEMENT << nl;

    os.writeKeyword("master") << master_ << token::END_STATEMENT << nl;

    if (dict_.found("note"))
    {
        os.writeKeyword("note") << string(dict_.lookup("note"))
            << token::END_STATEMENT << nl;
    }

    dict_.write(os);
}


bool Foam::interRegionHeatTransferModel::read(const dictionary& dict)
{
    if (basicSource::read(dict))
    {
        return true;
    }
    else
    {
        return false;
    }
}
// ************************************************************************* //
