/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

#include "constSolidThermo.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(constSolidThermo, 0);
    addToRunTimeSelectionTable(basicSolidThermo, constSolidThermo, mesh);
    addToRunTimeSelectionTable(basicSolidThermo, constSolidThermo, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constSolidThermo::constSolidThermo
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    basicSolidThermo(mesh, dict),
    dict_(dict.subDict(typeName + "Coeffs")),
    constKappa_(dimensionedScalar(dict_.lookup("kappa"))),
    kappa_
    (
        IOobject
        (
            "kappa",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        constKappa_
    ),
    constRho_(dimensionedScalar(dict_.lookup("rho"))),
    constCp_(dimensionedScalar(dict_.lookup("Cp"))),
    constHf_(dimensionedScalar(dict_.lookup("Hf"))),
    constEmissivity_(dimensionedScalar(dict_.lookup("emissivity"))),
    constKappaRad_(dimensionedScalar(dict_.lookup("kappaRad"))),
    constSigmaS_(dimensionedScalar(dict_.lookup("sigmaS")))
{
    read();

    kappa_ = constKappa_;

    rho_ = constRho_;

    emissivity_ = constEmissivity_;

    kappaRad_ = constKappaRad_;

    sigmaS_ = constSigmaS_;
}


Foam::constSolidThermo::constSolidThermo(const fvMesh& mesh)
:
    basicSolidThermo(mesh),
    dict_(subDict(typeName + "Coeffs")),
    constKappa_(dimensionedScalar(dict_.lookup("kappa"))),
    kappa_
    (
        IOobject
        (
            "kappa",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        constKappa_
    ),
    constRho_(dimensionedScalar(dict_.lookup("rho"))),
    constCp_(dimensionedScalar(dict_.lookup("Cp"))),
    constHf_(dimensionedScalar(dict_.lookup("Hf"))),
    constEmissivity_(dimensionedScalar(dict_.lookup("emissivity"))),
    constKappaRad_(dimensionedScalar(dict_.lookup("kappaRad"))),
    constSigmaS_(dimensionedScalar(dict_.lookup("sigmaS")))
{
    read();

    kappa_ = constKappa_;

    rho_ = constRho_;

    emissivity_ = constEmissivity_;

    kappaRad_ = constKappaRad_;

    sigmaS_ = constSigmaS_;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::constSolidThermo::~constSolidThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::constSolidThermo::correct()
{}


Foam::tmp<Foam::volScalarField> Foam::constSolidThermo::kappa() const
{
    return kappa_;
}


Foam::tmp<Foam::volSymmTensorField>
Foam::constSolidThermo::directionalKappa() const
{
    dimensionedSymmTensor t
    (
        constKappa_.name(),
        constKappa_.dimensions(),
        symmTensor
        (
            constKappa_.value(),
            0.0,
            0.0,
            constKappa_.value(),
            0.0,
            constKappa_.value()
        )
    );
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "kappa",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            t
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::constSolidThermo::Cp() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Cp",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            constCp_
        )
    );
}


Foam::tmp<Foam::volScalarField> Foam::constSolidThermo::Hf() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Hf",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            constHf_
        )
    );
}


Foam::tmp<Foam::scalarField> Foam::constSolidThermo::kappa
(
    const label patchI
) const
{
    return (kappa_.boundaryField()[patchI]);
}


Foam::tmp<Foam::symmTensorField> Foam::constSolidThermo::directionalKappa
(
    const label patchI
) const
{
    symmTensor t
    (
        constKappa_.value(),
        0.0,
        0.0,
        constKappa_.value(),
        0.0,
        constKappa_.value()
    );
    return tmp<symmTensorField>
    (
        new symmTensorField
        (
            T_.boundaryField()[patchI].size(),
            t
        )
    );
}


Foam::tmp<Foam::scalarField> Foam::constSolidThermo::Cp
(
    const label patchI
) const
{
    return tmp<scalarField>
    (
        new scalarField
        (
            T_.boundaryField()[patchI].size(),
            constCp_.value()
        )
    );
}


Foam::tmp<Foam::scalarField> Foam::constSolidThermo::Hf
(
    const label patchI
) const
{
    return tmp<scalarField>
    (
        new scalarField
        (
            T_.boundaryField()[patchI].size(),
            constHf_.value()
        )
    );
}


bool Foam::constSolidThermo::read()
{
    return read(subDict(typeName + "Coeffs"));
}


bool Foam::constSolidThermo::read(const dictionary& dict)
{
    constRho_ = dimensionedScalar(dict.lookup("rho"));
    constCp_ = dimensionedScalar(dict.lookup("Cp"));
    constKappa_ = dimensionedScalar(dict.lookup("kappa"));
    constHf_ = dimensionedScalar(dict.lookup("Hf"));
    constEmissivity_ = dimensionedScalar(dict.lookup("emissivity"));
    constKappaRad_ = dimensionedScalar(dict_.lookup("kappaRad"));
    constSigmaS_ = dimensionedScalar(dict_.lookup("sigmaS"));

    Info<< "Constructed constSolidThermo with" << nl
        << "    rho        : " << constRho_ << nl
        << "    Cp         : " << constCp_ << nl
        << "    kappa      : " << constKappa_ << nl
        << "    Hf         : " << constHf_ << nl
        << "    emissivity : " << constEmissivity_ << nl
        << "    kappaRad   : " << constKappaRad_ << nl
        << "    sigmaS     : " << constSigmaS_ << nl
        << endl;

    return true;
}


bool Foam::constSolidThermo::writeData(Ostream& os) const
{
    bool ok = basicSolidThermo::writeData(os);
    os.writeKeyword("rho") << constRho_ << token::END_STATEMENT << nl;
    os.writeKeyword("Cp") << constCp_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappa") << constKappa_ << token::END_STATEMENT << nl;
    os.writeKeyword("Hf") << constHf_ << token::END_STATEMENT << nl;
    os.writeKeyword("kappaRad") << constKappaRad_ << token::END_STATEMENT << nl;
    os.writeKeyword("sigmaS") << constSigmaS_ << token::END_STATEMENT << nl;
    os.writeKeyword("emissivity") << constEmissivity_ << token::END_STATEMENT
        << nl;
    return ok && os.good();
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const constSolidThermo& s)
{
    s.writeData(os);
    return os;
}


// ************************************************************************* //
