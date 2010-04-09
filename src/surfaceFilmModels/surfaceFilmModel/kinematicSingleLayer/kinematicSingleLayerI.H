/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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

#include "kinematicSingleLayer.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const Foam::volVectorField&
Foam::surfaceFilmModels::kinematicSingleLayer::nHat() const
{
    return nHat_;
}


inline const Foam::volScalarField&
Foam::surfaceFilmModels::kinematicSingleLayer::magSf() const
{
    return magSf_;
}


inline const Foam::Switch&
Foam::surfaceFilmModels::kinematicSingleLayer::momentumPredictor() const
{
    return momentumPredictor_;
}


inline Foam::label
Foam::surfaceFilmModels::kinematicSingleLayer::nOuterCorr() const
{
    return nOuterCorr_;
}


inline Foam::label Foam::surfaceFilmModels::kinematicSingleLayer::nCorr() const
{
    return nCorr_;
}


inline Foam::label
Foam::surfaceFilmModels::kinematicSingleLayer::nNonOrthCorr() const
{
    return nNonOrthCorr_;
}


inline Foam::scalar Foam::surfaceFilmModels::kinematicSingleLayer::Cf() const
{
    return Cf_;
}


inline const Foam::volScalarField&
Foam::surfaceFilmModels::kinematicSingleLayer::mu() const
{
    return mu_;
}


inline const Foam::volScalarField&
Foam::surfaceFilmModels::kinematicSingleLayer::sigma() const
{
    return sigma_;
}


inline const Foam::volScalarField&
Foam::surfaceFilmModels::kinematicSingleLayer::delta() const
{
    return delta_;
}

inline Foam::volVectorField&
Foam::surfaceFilmModels::kinematicSingleLayer::USpPrimary()
{
    return USpPrimary_;
}


inline Foam::volScalarField&
Foam::surfaceFilmModels::kinematicSingleLayer::pSpPrimary()
{
    return pSpPrimary_;
}


inline Foam::volScalarField&
Foam::surfaceFilmModels::kinematicSingleLayer::rhoSpPrimary()
{
    return rhoSpPrimary_;
}


inline Foam::volVectorField&
Foam::surfaceFilmModels::kinematicSingleLayer::USp()
{
    return USp_;
}


inline Foam::volScalarField&
Foam::surfaceFilmModels::kinematicSingleLayer::pSp()
{
    return pSp_;
}


inline Foam::volScalarField&
Foam::surfaceFilmModels::kinematicSingleLayer::rhoSp()
{
    return rhoSp_;
}


inline Foam::surfaceFilmModels::injectionModel&
Foam::surfaceFilmModels::kinematicSingleLayer::injection()
{
    return injection_();
}


inline Foam::tmp<Foam::volScalarField>
Foam::surfaceFilmModels::kinematicSingleLayer::gNorm() const
{
    tmp<volScalarField> tgNorm
    (
        new volScalarField
        (
            IOobject
            (
                "gNorm",
                filmRegion_.time().timeName(),
                filmRegion_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            g_ & nHat_,
            zeroGradientFvPatchScalarField::typeName
        )
    );

    return tgNorm;
}


inline Foam::tmp<Foam::volScalarField>
Foam::surfaceFilmModels::kinematicSingleLayer::gNormClipped() const
{
    tmp<volScalarField> tgNormClipped
    (
        new volScalarField
        (
            IOobject
            (
                "gNormClipped",
                filmRegion_.time().timeName(),
                filmRegion_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            g_ & nHat_,
            zeroGradientFvPatchScalarField::typeName
        )
    );

    volScalarField& gNormClipped = tgNormClipped();
    gNormClipped.min(0.0);

    return tgNormClipped;
}


inline Foam::tmp<Foam::volVectorField>
Foam::surfaceFilmModels::kinematicSingleLayer::gTan() const
{
    tmp<volVectorField> tgTan
    (
        new volVectorField
        (
            IOobject
            (
                "gTan",
                filmRegion_.time().timeName(),
                filmRegion_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            g_ - nHat_*gNorm(),
            zeroGradientFvPatchVectorField::typeName
        )
    );

    return tgTan;
}


// ************************************************************************* //
