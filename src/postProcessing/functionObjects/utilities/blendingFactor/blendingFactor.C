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

#include "blendingFactor.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(blendingFactor, 0);
}


// * * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * //

void Foam::blendingFactor::writeFileHeader(Ostream& os) const
{
    writeHeader(os, "Blending factor");
    writeCommented(os, "Time");
    writeTabbed(os, "Scheme1");
    writeTabbed(os, "Scheme2");
    writeTabbed(os, "Blended");
    os  << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blendingFactor::blendingFactor
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    functionObjectState(obr, name),
    functionObjectFile(obr, name, typeName, dict),
    name_(name),
    obr_(obr),
    phiName_("phi"),
    fieldName_("unknown-fieldName"),
    resultName_(word::null),
    tolerance_(0.001),
    log_(true)
{
    // Check if the available mesh is an fvMesh, otherwise deactivate
    if (setActive<fvMesh>())
    {
        read(dict);
        writeFileHeader(file());

        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        volScalarField* indicatorPtr
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
                dimensionedScalar("0", dimless, 0.0),
                zeroGradientFvPatchScalarField::typeName
            )
        );

        mesh.objectRegistry::store(indicatorPtr);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blendingFactor::~blendingFactor()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::blendingFactor::read(const dictionary& dict)
{
    if (active_)
    {
        functionObjectFile::read(dict);

        log_.readIfPresent("log", dict);

        dict.readIfPresent("phiName", phiName_);
        dict.lookup("fieldName") >> fieldName_;


        if (!dict.readIfPresent("resultName", resultName_))
        {
            resultName_ = name_ + ':' + fieldName_;
        }

        dict.readIfPresent("tolerance", tolerance_);
        if ((tolerance_ < 0) || (tolerance_ > 1))
        {
            FatalErrorInFunction
                << "tolerance must be in the range 0 to 1.  Supplied value: "
                << tolerance_ << exit(FatalError);
        }
    }
}


void Foam::blendingFactor::execute()
{
    if (active_)
    {
        calc<scalar>();
        calc<vector>();
    }
}


void Foam::blendingFactor::end()
{
    // Do nothing
}


void Foam::blendingFactor::timeSet()
{
    // Do nothing
}


void Foam::blendingFactor::write()
{
    if (active_)
    {
        const volScalarField& indicator =
            obr_.lookupObject<volScalarField>(resultName_);

        if (log_) Info
            << type() << " " << name_ << " output:" << nl
            << "    writing field " << indicator.name() << nl
            << endl;

        indicator.write();
    }
}


// ************************************************************************* //
