/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd.
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

#include "DESModelRegions.H"
#include "volFields.H"
#include "DESModelBase.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(DESModelRegions, 0);
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::DESModelRegions::writeFileHeader(Ostream& os) const
{
    writeHeader(os, "DES model region coverage (% volume)");

    writeCommented(os, "Time");
    writeTabbed(os, "LES");
    writeTabbed(os, "RAS");
    os << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::DESModelRegions::DESModelRegions
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
    resultName_(name),
    log_(true)
{
    // Check if the available mesh is an fvMesh, otherwise deactivate
    if (!isA<fvMesh>(obr_))
    {
        active_ = false;
        WarningIn
        (
            "DESModelRegions::DESModelRegions"
            "("
                "const word&, "
                "const objectRegistry&, "
                "const dictionary&, "
                "const bool"
            ")"
        )   << "No fvMesh available, deactivating " << name_ << nl
            << endl;
    }

    read(dict);

    if (active_)
    {
        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        volScalarField* DESModelRegionsPtr
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

        mesh.objectRegistry::store(DESModelRegionsPtr);

        writeFileHeader(file());
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::DESModelRegions::~DESModelRegions()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::DESModelRegions::read(const dictionary& dict)
{
    if (active_)
    {
        functionObjectFile::read(dict);

        log_.readIfPresent("log", dict);
        dict.readIfPresent("resultName", resultName_);
    }
}


void Foam::DESModelRegions::execute()
{
    if (active_)
    {
        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        if (log_) Info<< type() << " " << name_ <<  " output:" << nl;

        volScalarField& DESModelRegions =
            const_cast<volScalarField&>
            (
                mesh.lookupObject<volScalarField>(resultName_)
            );


        if (mesh.foundObject<DESModelBase>(turbulenceModel::propertiesName))
        {
            const DESModelBase& model =
                mesh.lookupObject<DESModelBase>
                (
                    turbulenceModel::propertiesName
                );

            DESModelRegions == model.LESRegion();

            scalar prc =
                gSum(DESModelRegions.internalField()*mesh.V())
               /gSum(mesh.V())*100.0;

            file() << obr_.time().value()
                << token::TAB << prc
                << token::TAB << 100.0 - prc
                << endl;

            if (log_) Info
                << "    LES = " << prc << " % (volume)" << nl
                << "    RAS = " << 100.0 - prc << " % (volume)" << nl
                << endl;
        }
        else
        {
            if (log_) Info
                << "    No DES turbulence model found in database" << nl
                << endl;
        }
    }
}


void Foam::DESModelRegions::end()
{
    // Do nothing
}


void Foam::DESModelRegions::timeSet()
{
    // Do nothing
}


void Foam::DESModelRegions::write()
{
    if (active_)
    {
        const volScalarField& DESModelRegions =
            obr_.lookupObject<volScalarField>(resultName_);

        if (log_) Info
            << type() << " " << name_ <<  " output:" << nl
            << "    writing field " << DESModelRegions.name() << nl
            << endl;

        DESModelRegions.write();
    }
}


// ************************************************************************* //
