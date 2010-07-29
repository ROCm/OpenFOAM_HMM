/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

#include "interpolatedSolidThermo.H"
#include "addToRunTimeSelectionTable.H"
#include "interpolateXY.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(interpolatedSolidThermo, 0);
    addToRunTimeSelectionTable
    (
        basicSolidThermo,
        interpolatedSolidThermo,
        mesh
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interpolatedSolidThermo::interpolatedSolidThermo(const fvMesh& mesh)
:
    basicSolidThermo(mesh)
{
    read();
    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interpolatedSolidThermo::~interpolatedSolidThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::interpolatedSolidThermo::correct()
{}


Foam::tmp<Foam::volScalarField> Foam::interpolatedSolidThermo::rho() const
{
    tmp<volScalarField> trho
    (
        new volScalarField
        (
            IOobject
            (
                "rho",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimDensity
        )
    );
    volScalarField& rho = trho();

    rho.internalField() = interpolateXY
    (
        T_.internalField(),
        TValues_,
        rhoValues_
    );

    forAll(rho.boundaryField(), patchI)
    {
        rho.boundaryField()[patchI] == this->rho(patchI)();
    }

    return trho;
}


Foam::tmp<Foam::volScalarField> Foam::interpolatedSolidThermo::cp() const
{
    tmp<volScalarField> tcp
    (
        new volScalarField
        (
            IOobject
            (
                "cp",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimEnergy/(dimMass*dimTemperature)
        )
    );
    volScalarField& cp = tcp();

    cp.internalField() = interpolateXY
    (
        T_.internalField(),
        TValues_,
        cpValues_
    );

    forAll(cp.boundaryField(), patchI)
    {
        cp.boundaryField()[patchI] == this->cp(patchI)();
    }

    return tcp;
}


Foam::tmp<Foam::volScalarField> Foam::interpolatedSolidThermo::K() const
{
    tmp<volScalarField> tK
    (
        new volScalarField
        (
            IOobject
            (
                "K",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimEnergy/dimTime/(dimLength*dimTemperature)
        )
    );
    volScalarField& K = tK();

    K.internalField() = interpolateXY
    (
        T_.internalField(),
        TValues_,
        KValues_
    );

    forAll(K.boundaryField(), patchI)
    {
        K.boundaryField()[patchI] == this->K(patchI)();
    }

    return tK;
}


Foam::tmp<Foam::volSymmTensorField>
Foam::interpolatedSolidThermo::directionalK()
const
{
    tmp<volSymmTensorField> tK
    (
        new volSymmTensorField
        (
            IOobject
            (
                "K",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedSymmTensor
            (
                "zero",
                dimEnergy/dimTime/(dimLength*dimTemperature),
                symmTensor::zero
            )
        )
    );
    volSymmTensorField& K = tK();

    Field<scalar> scalarK
    (
        interpolateXY
        (
            T_.internalField(),
            TValues_,
            KValues_
        )
    );

    K.internalField().replace(symmTensor::XX, scalarK);
    K.internalField().replace(symmTensor::YY, scalarK);
    K.internalField().replace(symmTensor::ZZ, scalarK);

    forAll(K.boundaryField(), patchI)
    {
        K.boundaryField()[patchI] == this->directionalK(patchI)();
    }

    return tK;
}


Foam::tmp<Foam::volScalarField> Foam::interpolatedSolidThermo::Hf() const
{
    tmp<volScalarField> tHf
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
            dimEnergy/dimMass
        )
    );
    volScalarField& Hf = tHf();

    Hf.internalField() = interpolateXY
    (
        T_.internalField(),
        TValues_,
        HfValues_
    );

    forAll(Hf.boundaryField(), patchI)
    {
        Hf.boundaryField()[patchI] == this->Hf(patchI)();
    }

    return tHf;
}


Foam::tmp<Foam::volScalarField>
Foam::interpolatedSolidThermo::emissivity() const
{
    tmp<volScalarField> temissivity
    (
        new volScalarField
        (
            IOobject
            (
                "emissivity",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimless
        )
    );
    volScalarField& emissivity = temissivity();

    emissivity.internalField() = interpolateXY
    (
        T_.internalField(),
        TValues_,
        emissivityValues_
    );

    forAll(emissivity.boundaryField(), patchI)
    {
        emissivity.boundaryField()[patchI] == this->emissivity(patchI)();
    }

    return temissivity;
}


Foam::tmp<Foam::scalarField> Foam::interpolatedSolidThermo::rho
(
    const label patchI
) const
{
    return tmp<scalarField>
    (
        new scalarField
        (
            interpolateXY
            (
                T_.boundaryField()[patchI],
                TValues_,
                rhoValues_
            )
        )
    );
}


Foam::tmp<Foam::scalarField> Foam::interpolatedSolidThermo::cp
(
    const label patchI
) const
{
    return tmp<scalarField>
    (
        new scalarField
        (
            interpolateXY
            (
                T_.boundaryField()[patchI],
                TValues_,
                cpValues_
            )
        )
    );
}


Foam::tmp<Foam::scalarField> Foam::interpolatedSolidThermo::K
(
    const label patchI
) const
{
    return tmp<scalarField>
    (
        new scalarField
        (
            interpolateXY
            (
                T_.boundaryField()[patchI],
                TValues_,
                KValues_
            )
        )
    );
}


Foam::tmp<Foam::symmTensorField> Foam::interpolatedSolidThermo::directionalK
(
    const label patchI
) const
{
    const fvPatchScalarField& patchT = T_.boundaryField()[patchI];

    Field<scalar> scalarK(interpolateXY(patchT, TValues_, KValues_));

    tmp<symmTensorField> tfld
    (
        new symmTensorField
        (
            scalarK.size(),
            symmTensor::zero
        )
    );
    symmTensorField& fld = tfld();

    fld.replace(symmTensor::XX, scalarK);
    fld.replace(symmTensor::YY, scalarK);
    fld.replace(symmTensor::ZZ, scalarK);

    return tfld;
}


Foam::tmp<Foam::scalarField> Foam::interpolatedSolidThermo::Hf
(
    const label patchI
) const
{
    return tmp<scalarField>
    (
        new scalarField
        (
            interpolateXY
            (
                T_.boundaryField()[patchI],
                TValues_,
                HfValues_
            )
        )
    );
}


Foam::tmp<Foam::scalarField> Foam::interpolatedSolidThermo::emissivity
(
    const label patchI
) const
{
    return tmp<scalarField>
    (
        new scalarField
        (
            interpolateXY
            (
                T_.boundaryField()[patchI],
                TValues_,
                emissivityValues_
            )
        )
    );
}


bool Foam::interpolatedSolidThermo::read()
{
    return read(subDict(typeName + "Coeffs"));
}


bool Foam::interpolatedSolidThermo::read(const dictionary& dict)
{
    TValues_ = Field<scalar>(dict.lookup("TValues"));
    rhoValues_ = Field<scalar>(dict.lookup("rhoValues"));
    cpValues_ = Field<scalar>(dict.lookup("cpValues"));
    KValues_ = Field<scalar>(dict.lookup("KValues"));
    HfValues_ = Field<scalar>(dict.lookup("HfValues"));
    emissivityValues_ = Field<scalar>(dict.lookup("emissivityValues"));

    Info<< "Constructed interpolatedSolidThermo with samples" << nl
        << "    T          : " << TValues_ << nl
        << "    rho        : " << rhoValues_ << nl
        << "    cp         : " << cpValues_ << nl
        << "    K          : " << KValues_ << nl
        << "    Hf         : " << HfValues_ << nl
        << "    emissivity : " << emissivityValues_ << nl
        << endl;

    if
    (
        (TValues_.size() != rhoValues_.size())
     && (TValues_.size() != cpValues_.size())
     && (TValues_.size() != rhoValues_.size())
     && (TValues_.size() != KValues_.size())
     && (TValues_.size() != HfValues_.size())
     && (TValues_.size() != emissivityValues_.size())
    )
    {
        FatalIOErrorIn("interpolatedSolidThermo::read()", dict)
            << "Size of property tables should be equal to size of Temperature"
            << " values " << TValues_.size()
            << exit(FatalIOError);
    }

    for (label i = 1; i < TValues_.size(); i++)
    {
        if (TValues_[i] <= TValues_[i-1])
        {
            FatalIOErrorIn("interpolatedSolidThermo::read()", dict)
                << "Temperature values are not in increasing order "
                << TValues_ << exit(FatalIOError);
        }
    }
    return true;
}


bool Foam::interpolatedSolidThermo::writeData(Ostream& os) const
{
    bool ok = basicSolidThermo::writeData(os);
    os.writeKeyword("TValues") << TValues_ << token::END_STATEMENT << nl;
    os.writeKeyword("rhoValues") << rhoValues_ << token::END_STATEMENT << nl;
    os.writeKeyword("cpValues") << cpValues_ << token::END_STATEMENT << nl;
    os.writeKeyword("KValues") << KValues_ << token::END_STATEMENT << nl;
    os.writeKeyword("HfValues") << HfValues_ << token::END_STATEMENT << nl;
    os.writeKeyword("emissivityValues") << emissivityValues_
        << token::END_STATEMENT << nl;

    return ok && os.good();
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const interpolatedSolidThermo& s)
{
    s.writeData(os);
    return os;
}


// ************************************************************************* //
