/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "filmThermoModel.H"
#include "surfaceInterpolate.H"
#include "fvcSurfaceIntegrate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace surfaceFilmModels
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline Switch kinematicSingleLayer::momentumPredictor() const
{
    return momentumPredictor_;
}


inline label kinematicSingleLayer::nOuterCorr() const
{
    return nOuterCorr_;
}


inline label kinematicSingleLayer::nCorr() const
{
    return nCorr_;
}


inline label kinematicSingleLayer::nNonOrthCorr() const
{
    return nNonOrthCorr_;
}


inline const dimensionedScalar& kinematicSingleLayer::deltaSmall() const
{
    return deltaSmall_;
}


inline const volScalarField& kinematicSingleLayer::mu() const
{
    return mu_;
}


inline const volScalarField& kinematicSingleLayer::sigma() const
{
    return sigma_;
}


inline const volScalarField& kinematicSingleLayer::delta() const
{
    return delta_;
}


inline const volScalarField& kinematicSingleLayer::alpha() const
{
    return alpha_;
}


inline volVectorField& kinematicSingleLayer::USpPrimary()
{
    return USpPrimary_;
}


inline volScalarField& kinematicSingleLayer::pSpPrimary()
{
    return pSpPrimary_;
}


inline volScalarField& kinematicSingleLayer::rhoSpPrimary()
{
    return rhoSpPrimary_;
}


inline volVectorField& kinematicSingleLayer::USp()
{
    return USp_;
}


inline volScalarField& kinematicSingleLayer::pSp()
{
    return pSp_;
}


inline volScalarField& kinematicSingleLayer::rhoSp()
{
    return rhoSp_;
}


inline const volVectorField& kinematicSingleLayer::USp() const
{
    return USp_;
}


inline const volScalarField& kinematicSingleLayer::pSp() const
{
    return pSp_;
}


inline const volScalarField& kinematicSingleLayer::rhoSp() const
{
    return rhoSp_;
}


inline const volVectorField& kinematicSingleLayer::UPrimary() const
{
    return UPrimary_;
}


inline const volScalarField& kinematicSingleLayer::pPrimary() const
{
    return pPrimary_;
}


inline const volScalarField& kinematicSingleLayer::rhoPrimary() const
{
    return rhoPrimary_;
}


inline const volScalarField& kinematicSingleLayer::muPrimary() const
{
    return muPrimary_;
}


inline const filmThermoModel& kinematicSingleLayer::filmThermo() const
{
    return *filmThermo_;
}


inline injectionModelList& kinematicSingleLayer::injection()
{
    return injection_;
}


inline transferModelList& kinematicSingleLayer::transfer()
{
    return transfer_;
}


inline const filmTurbulenceModel& kinematicSingleLayer::turbulence() const
{
    return *turbulence_;
}


inline tmp<volScalarField> kinematicSingleLayer::mass() const
{
    return deltaRho_*magSf();
}


inline tmp<volScalarField> kinematicSingleLayer::deltaMass() const
{
    return rhoSp_*magSf()*time().deltaT();
}


inline tmp<volScalarField> kinematicSingleLayer::gNorm() const
{
    tmp<volScalarField> tgNorm
    (
        new volScalarField
        (
            IOobject
            (
                "gNorm",
                time().timeName(),
                regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            g_ & nHat()
        )
    );

    return tgNorm;
}


inline tmp<volScalarField> kinematicSingleLayer::gNormClipped() const
{
    tmp<volScalarField> tgNormClipped
    (
        new volScalarField
        (
            IOobject
            (
                "gNormClipped",
                time().timeName(),
                regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            g_ & nHat()
        )
    );

    volScalarField& gNormClipped = tgNormClipped.ref();
    gNormClipped.min(0.0);

    return tgNormClipped;
}


inline tmp<volVectorField> kinematicSingleLayer::gTan() const
{
    tmp<volVectorField> tgTan
    (
        new volVectorField
        (
            IOobject
            (
                "gTan",
                time().timeName(),
                regionMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            g_ - nHat()*gNorm()
        )
    );

    return tgTan;
}

inline tmp<vectorField> kinematicSingleLayer::gTan
(
    const label patchI
) const
{
    const vectorField& nH = nHat().boundaryField()[patchI];
    const vector& g = g_.value();
    tmp<vectorField> tgTan(new vectorField(g - nH*(g & nH)));

    return tgTan;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace surfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
