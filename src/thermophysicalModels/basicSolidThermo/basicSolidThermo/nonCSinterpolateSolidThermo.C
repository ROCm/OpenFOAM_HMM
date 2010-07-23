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

#include "interpolateSolidThermo.H"
#include "addToRunTimeSelectionTable.H"
#include "interpolateXY.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(interpolateSolidThermo, 0);
    addToRunTimeSelectionTable
    (
        basicSolidThermo,
        interpolateSolidThermo,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interpolateSolidThermo::interpolateSolidThermo
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    basicSolidThermo(mesh, dict, typeName),
    TValues_(dict_.lookup("TValues")),
    rhoValues_(dict_.lookup("rhoValues")),
    cpValues_(dict_.lookup("cpValues")),
    KValues_(dict_.lookup("KValues")),
    HfValues_(dict_.lookup("HfValues")),
    emissivityValues_(dict_.lookup("emissivityValues"))
{
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
        FatalIOErrorIn
        (
            "interpolateSolidThermo::interpolateSolidThermo\n"
            "(\n"
            "    const fvMesh& mesh,\n"
            "    const dictionary& dict\n"
            ")\n",
            dict_
        )   << "Size of property tables should be equal to size of Temperature"
            << " values " << TValues_.size()
            << exit(FatalIOError);
    }

    for (label i = 1; i < TValues_.size(); i++)
    {
        if (TValues_[i] <= TValues_[i-1])
        {
            FatalIOErrorIn
            (
                "interpolateSolidThermo::interpolateSolidThermo\n"
                "(\n"
                "    const fvMesh& mesh,\n"
                "    const dictionary& dict\n"
                ")\n",
                dict_
            )   << "Temperature values are not in increasing order "
                << TValues_ << exit(FatalIOError);
        }
    }

    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::interpolateSolidThermo::~interpolateSolidThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::interpolateSolidThermo::correct()
{
    // rho
    rho_.internalField() = interpolateXY
    (
        T_.internalField(),
        TValues_,
        rhoValues_
    );

    forAll(rho_.boundaryField(), patchI)
    {
        rho_.boundaryField()[patchI] == interpolateXY
        (
            T_.boundaryField()[patchI],
            TValues_,
            rhoValues_
        );
    }


    // cp
    cp_.internalField() = interpolateXY
    (
        T_.internalField(),
        TValues_,
        cpValues_
    );

    forAll(cp_.boundaryField(), patchI)
    {
        cp_.boundaryField()[patchI] == interpolateXY
        (
            T_.boundaryField()[patchI],
            TValues_,
            cpValues_
        );
    }


    // K
    K_.internalField() = interpolateXY
    (
        T_.internalField(),
        TValues_,
        KValues_
    );

    forAll(K_.boundaryField(), patchI)
    {
        K_.boundaryField()[patchI] == interpolateXY
        (
            T_.boundaryField()[patchI],
            TValues_,
            KValues_
        );
    }


    // Hf
    Hf_.internalField() = interpolateXY
    (
        T_.internalField(),
        TValues_,
        HfValues_
    );

    forAll(Hf_.boundaryField(), patchI)
    {
        Hf_.boundaryField()[patchI] == interpolateXY
        (
            T_.boundaryField()[patchI],
            TValues_,
            HfValues_
        );
    }


    // emissivity
    emissivity_.internalField() = interpolateXY
    (
        T_.internalField(),
        TValues_,
        emissivityValues_
    );

    forAll(emissivity_.boundaryField(), patchI)
    {
        emissivity_.boundaryField()[patchI] == interpolateXY
        (
            T_.boundaryField()[patchI],
            TValues_,
            emissivityValues_
        );
    }
}


void Foam::interpolateSolidThermo::write(Ostream& os) const
{
    basicSolidThermo::write(os);
    os.writeKeyword("TValues") << TValues_ << token::END_STATEMENT << nl;
    os.writeKeyword("rhoValues") << rhoValues_ << token::END_STATEMENT << nl;
    os.writeKeyword("cpValues") << cpValues_ << token::END_STATEMENT << nl;
    os.writeKeyword("KValues") << KValues_ << token::END_STATEMENT << nl;
    os.writeKeyword("HfValues") << HfValues_ << token::END_STATEMENT << nl;
    os.writeKeyword("emissivityValues") << emissivityValues_
        << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const interpolateSolidThermo& s)
{
    s.write(os);
    return os;
}


// ************************************************************************* //
