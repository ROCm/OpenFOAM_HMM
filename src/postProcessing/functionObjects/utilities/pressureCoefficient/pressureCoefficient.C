/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

#include "pressureCoefficient.H"
#include "volFields.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::pressureCoefficient, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::pressureCoefficient::rho() const
{
    if (rhoName_ == "rhoInf")
    {
        const fvMesh& mesh = refCast<const fvMesh>(obr_);

        return tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "rho",
                    mesh.time().timeName(),
                    mesh
                ),
                mesh,
                dimensionedScalar("rho", dimDensity, rhoRef_)
            )
        );
    }
    else
    {
        return(obr_.lookupObject<volScalarField>(rhoName_));
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pressureCoefficient::pressureCoefficient
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    active_(true),
    pName_("p"),
    rhoName_(word::null),
    rhoRef_(-GREAT),
    magUinf_(0.0)
{
    // Check if the available mesh is an fvMesh, otherwise deactivate
    if (!isA<fvMesh>(obr_))
    {
        active_ = false;
        WarningIn
        (
            "pressureCoefficient::pressureCoefficient"
            "("
                "const word&, "
                "const objectRegistry&, "
                "const dictionary&, "
                "const bool"
            ")"
        )   << "No fvMesh available, deactivating." << nl
            << endl;
    }

    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pressureCoefficient::~pressureCoefficient()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pressureCoefficient::read(const dictionary& dict)
{
    if (active_)
    {
        pName_ = dict.lookupOrDefault<word>("pName", "p");

        rhoName_ = dict.lookupOrDefault<word>("rhoName", "rho");
        if (rhoName_ == "rhoInf")
        {
            dict.lookup("rhoInf") >> rhoRef_;
        }

        dict.lookup("magUinf") >> magUinf_;
    }
}


void Foam::pressureCoefficient::execute()
{
    // Do nothing - only valid on write
}


void Foam::pressureCoefficient::end()
{
    // Do nothing - only valid on write
}


void Foam::pressureCoefficient::write()
{
    if (active_)
    {
        const volScalarField& p = obr_.lookupObject<volScalarField>(pName_);

        volScalarField pressureCoefficient
        (
            IOobject
            (
                "pressureCoefficient",
                obr_.time().timeName(),
                obr_,
                IOobject::NO_READ
            ),
            p/(0.5*rho()*sqr(magUinf_))
        );

        pressureCoefficient.write();
    }
}


// ************************************************************************* //
