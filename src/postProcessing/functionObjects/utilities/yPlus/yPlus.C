/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd
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

#include "yPlus.H"
#include "volFields.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(yPlus, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::yPlus::writeFileHeader(Ostream& os) const
{
    writeHeader(os, "y+");
    writeCommented(os, "Time");
    writeTabbed(os, "patch");
    writeTabbed(os, "min");
    writeTabbed(os, "max");
    writeTabbed(os, "average");
    os << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::yPlus::yPlus
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    functionObjectFile(obr, name, typeName, dict),
    name_(name),
    obr_(obr),
    active_(true),
    phiName_("phi"),
    resultName_(name),
    log_(true)
{
    // Check if the available mesh is an fvMesh, otherwise deactivate
    if (!isA<fvMesh>(obr_))
    {
        active_ = false;
        WarningIn
        (
            "yPlus::yPlus"
            "("
                "const word&, "
                "const objectRegistry&, "
                "const dictionary&, "
                "const bool"
            ")"
        )   << "No fvMesh available, deactivating " << name_ << nl
            << endl;
    }

    if (active_)
    {
        read(dict);

        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        volScalarField* yPlusPtr
        (
            new volScalarField
            (
                IOobject
                (
                    resultName_,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh,
                dimensionedScalar("0", dimless, 0.0)
            )
        );

        mesh.objectRegistry::store(yPlusPtr);

        writeFileHeader(file());
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::yPlus::~yPlus()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::yPlus::read(const dictionary& dict)
{
    if (active_)
    {
        functionObjectFile::read(dict);

        log_.readIfPresent("log", dict);
        dict.readIfPresent("resultName", resultName_);
        dict.readIfPresent("phiName", phiName_);
    }
}


void Foam::yPlus::execute()
{
    typedef compressible::turbulenceModel cmpTurbModel;
    typedef incompressible::turbulenceModel icoTurbModel;

    if (active_)
    {
        const surfaceScalarField& phi =
            obr_.lookupObject<surfaceScalarField>(phiName_);

        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        volScalarField& yPlus =
            const_cast<volScalarField&>
            (
                mesh.lookupObject<volScalarField>(resultName_)
            );

        if (log_) Info << type() << " " << name_ << " output:" << nl;

        if (phi.dimensions() == dimMass/dimTime)
        {
            if (mesh.foundObject<cmpTurbModel>(turbulenceModel::propertiesName))
            {
                const cmpTurbModel& model =
                    mesh.lookupObject<cmpTurbModel>
                    (
                        turbulenceModel::propertiesName
                    );

                calcYPlus(model, mesh, yPlus);
            }
            else
            {
                WarningIn("void Foam::yPlus::execute()")
                    << "Unable to find compressible turbulence model in the "
                    << "database: yPlus will not be calculated" << endl;
            }
        }
        else if (phi.dimensions() == dimVolume/dimTime)
        {
            if (mesh.foundObject<icoTurbModel>(turbulenceModel::propertiesName))
            {
                const icoTurbModel& model =
                    mesh.lookupObject<icoTurbModel>
                    (
                        turbulenceModel::propertiesName
                    );

                calcYPlus(model, mesh, yPlus);
            }
            else
            {
                WarningIn("void Foam::yPlus::execute()")
                    << "Unable to find incompressible turbulence model in the "
                    << "database: yPlus will not be calculated" << endl;
            }
        }
        else
        {
            WarningIn("void Foam::yPlus::execute()")
                << "Unknown " << phiName_ << " dimensions: "
                << phi.dimensions() << nl
                << "Expected either " << dimMass/dimTime << " or "
                << dimVolume/dimTime << nl
                << "yPlus will not be calculated" << endl;
        }
    }
}


void Foam::yPlus::end()
{
    // Do nothing
}


void Foam::yPlus::timeSet()
{
    // Do nothing
}


void Foam::yPlus::write()
{
    if (active_)
    {
        const volScalarField& yPlus =
            obr_.lookupObject<volScalarField>(resultName_);

        if (log_) Info
            << type() << " " << name_ << " output:" << nl
            << "    writing field " << yPlus.name() << nl
            << endl;

        yPlus.write();
    }
}


// ************************************************************************* //
