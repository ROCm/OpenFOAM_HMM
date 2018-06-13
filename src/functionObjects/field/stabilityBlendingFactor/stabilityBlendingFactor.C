/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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

#include "stabilityBlendingFactor.H"
#include "addToRunTimeSelectionTable.H"
#include "zeroGradientFvPatchFields.H"
#include "fvcGrad.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(stabilityBlendingFactor, 0);
    addToRunTimeSelectionTable
    (
        functionObject,
        stabilityBlendingFactor,
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::stabilityBlendingFactor::writeFileHeader
(
    Ostream& os
) const
{
    writeHeader(os, "Stabilization blending factor");
    writeCommented(os, "Time");
    writeTabbed(os, "Scheme1");
    writeTabbed(os, "Scheme2");
    writeTabbed(os, "Blended");
    os  << endl;
}

bool Foam::functionObjects::stabilityBlendingFactor::calc()
{
    init(true);
    return true;
}


bool Foam::functionObjects::stabilityBlendingFactor::init(bool first)
{
    const IOField<scalar>* residualPtr =
        mesh_.lookupObjectPtr<IOField<scalar>>(residualName_);

    if (residuals_)
    {
        if (!residualPtr)
        {
             WarningInFunction
                << "Could not find residual file : " << residualName_ << nl
                << "The residual mode won't be considered for the blended "
                << "field in the stability blending factor. " << nl
                << "Add the corresponding residual function object. " << nl
                << "If the residual function object is already set " << nl
                << "you might need to wait for the first iteration."
                << endl;
        }
        else
        {
            scalar meanRes = gAverage(mag(*residualPtr)) + VSMALL;

            if (log)
            {
                Log << "    Average(mag(residuals)) :  " << meanRes << endl;
            }

            oldError_ = error_;
            oldErrorIntegral_ = errorIntegral_;
            error_ = mag(meanRes - mag(*residualPtr));
            errorIntegral_ = oldErrorIntegral_ + 0.5*(error_ + oldError_);
            const scalarField errorDifferential(error_ - oldError_);

            const scalarField factorList
            (
                + P_*error_
                + I_*errorIntegral_
                + D_*errorDifferential
            );

            const scalarField indicatorResidual
            (
                max
                (
                    min
                    (
                        mag(factorList - meanRes)/(maxResidual_*meanRes),
                        1.0
                    ),
                    0.0
                )
            );

            forAll (indicator_, i)
            {
                indicator_[i] = indicatorResidual[i];
            }
        }
    }

    const volScalarField* nonOrthPtr =
        mesh_.lookupObjectPtr<volScalarField>(nonOrthogonalityName_);

    if (nonOrthogonality_)
    {
        indicator_ =
        max
        (
            indicator_,
            min
            (
                max
                (
                    0.0,
                    (*nonOrthPtr - maxNonOrthogonality_)
                   /(minNonOrthogonality_ - maxNonOrthogonality_)
                ),
                1.0
            )
        );

        if (log)
        {
            Log << "    Max non-orthogonality :  " << max(*nonOrthPtr).value()
                << endl;
        }
    }

    const volScalarField* skewnessPtr =
        mesh_.lookupObjectPtr<volScalarField>(skewnessName_);

    if (skewness_)
    {
        indicator_ =
        max
        (
            indicator_,
            min
            (
                max
                (
                    0.0,
                    (*skewnessPtr - maxSkewness_)
                  / (minSkewness_ - maxSkewness_)
                ),
                1.0
            )
        );

        if (log)
        {
            Log << "    Max skewness :  " << max(*skewnessPtr).value()
                << endl;
        }
    }

    const volScalarField* faceWeightsPtr =
        mesh_.lookupObjectPtr<volScalarField>(faceWeightName_);

    if (faceWeight_)
    {
        indicator_ =
            max
            (
                indicator_,
                min
                (
                    max
                    (
                        0.0,
                        (minFaceWeight_ - *faceWeightsPtr)
                      / (minFaceWeight_ - maxFaceWeight_)
                    ),
                    1.0
                )
            );

        if (log)
        {
            Log << "    Min face weight:  " << min(*faceWeightsPtr).value()
                << endl;
        }
    }


    if (gradCc_)
    {
        tmp<volScalarField> magGradCCPtr
        (
            new volScalarField
            (
                IOobject
                (
                    "magGradCC",
                    time_.timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar(dimless, Zero),
                zeroGradientFvPatchScalarField::typeName
            )
        );

        for (direction i=0; i<vector::nComponents; i++)
        {
            // Create field with zero grad
            volScalarField cci
            (
                IOobject
                (
                    "cc" + word(vector::componentNames[i]),
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar(dimLength, Zero),
                zeroGradientFvPatchScalarField::typeName
            );
            cci = mesh_.C().component(i);
            cci.correctBoundaryConditions();
            magGradCCPtr.ref() +=  mag(fvc::grad(cci)).ref();
        }

        if (log)
        {
            Log << "    Max magGradCc :  " << max(magGradCCPtr.ref()).value()
                << endl;
        }

        indicator_ =
            max
            (
                indicator_,
                min
                (
                    max
                    (
                        0.0,
                        (magGradCCPtr.ref() - maxGradCc_)
                      / (minGradCc_ - maxGradCc_)
                    ),
                    1.0
                )
            );
    }


    const volVectorField* UNamePtr =
        mesh_.lookupObjectPtr<volVectorField>(UName_);

    if (Co_)
    {
        tmp<volScalarField> CoPtr
        (
            new volScalarField
            (
                IOobject
                (
                    "Co",
                    time_.timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar(dimless, Zero),
                zeroGradientFvPatchScalarField::typeName
            )
        );

        volScalarField& Co = CoPtr.ref();

        Co.primitiveFieldRef() =
            mesh_.time().deltaT()*mag(*UNamePtr)/pow(mesh_.V(), 1.0/3.0);

        indicator_ =
            max
            (
                indicator_,
                min(max(0.0, (Co - Co1_)/(Co2_ - Co1_)), 1.0)
            );

        if (log)
        {
            Log << "    Max Co :  " << max(Co).value()
                << endl;
        }
    }

    indicator_.correctBoundaryConditions();
    indicator_.min(1.0);
    indicator_.max(0.0);

    // Update the blended surface field
    surfaceScalarField* surBlendedPtr =
    (
        mesh_.lookupObjectRefPtr<surfaceScalarField>(resultName_)
    );

    *surBlendedPtr = fvc::interpolate(indicator_);

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::stabilityBlendingFactor::stabilityBlendingFactor
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict),
    writeFile(obr_, name, typeName, dict),
    indicator_
    (
        IOobject
        (
            "blendedIndicator",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimless, Zero),
        zeroGradientFvPatchScalarField::typeName
    ),
    nonOrthogonality_(dict.lookupOrDefault<Switch>("switchNonOrtho", false)),
    gradCc_(dict.lookupOrDefault<Switch>("switchGradCc", false)),
    residuals_(dict.lookupOrDefault<Switch>("switchResiduals", false)),
    faceWeight_(dict.lookupOrDefault<Switch>("switchFaceWeight", false)),
    skewness_(dict.lookupOrDefault<Switch>("switchSkewness", false)),
    Co_(dict.lookupOrDefault<Switch>("switchCo", false)),

    maxNonOrthogonality_
    (
        dict.lookupOrDefault<scalar>("maxNonOrthogonality", 20.0)
    ),
    minNonOrthogonality_
    (
        dict.lookupOrDefault<scalar>("minNonOrthogonality", 60.0)
    ),
    maxGradCc_(dict.lookupOrDefault<scalar>("maxGradCc", 3.0)),
    minGradCc_(dict.lookupOrDefault<scalar>("minGradCc", 4.0)),
    maxResidual_(dict.lookupOrDefault<scalar>("maxResidual", 10.0)),
    minFaceWeight_(dict.lookupOrDefault<scalar>("minFaceWeight", 0.3)),
    maxFaceWeight_(dict.lookupOrDefault<scalar>("maxFaceWeight", 0.2)),
    maxSkewness_(dict.lookupOrDefault<scalar>("maxSkewness", 2.0)),
    minSkewness_(dict.lookupOrDefault<scalar>("minSkewness", 3.0)),
    Co1_(dict.lookupOrDefault<scalar>("Co1", 1.0)),
    Co2_(dict.lookupOrDefault<scalar>("Co2", 10.0)),

    nonOrthogonalityName_
    (
        dict.lookupOrDefault<word>("nonOrthogonality", "nonOrthoAngle")
    ),
    faceWeightName_
    (
        dict.lookupOrDefault<word>("faceWeight", "faceWeight")
    ),
    skewnessName_
    (
        dict.lookupOrDefault<word>("skewness", "skewness")
    ),
    residualName_
    (
        dict.lookupOrDefault<word>("residual", "initialResidual:p")
    ),
    UName_
    (
         dict.lookupOrDefault<word>("U", "U")
    ),

    tolerance_(0.001),
    error_(mesh_.nCells(), 0.0),
    errorIntegral_(mesh_.nCells(), 0.0),
    oldError_(mesh_.nCells(), 0.0),
    oldErrorIntegral_(mesh_.nCells(), 0.0),
    P_(dict.lookupOrDefault<scalar>("P", 3)),
    I_(dict.lookupOrDefault<scalar>("I", 0.0)),
    D_(dict.lookupOrDefault<scalar>("D", 0.25))
{
    read(dict);
    setResultName(typeName, "");

    tmp<surfaceScalarField> faceBlendedPtr
    (
        new surfaceScalarField
        (
            IOobject
            (
                resultName_,
                time_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar(dimless, Zero)
        )
    );
    store(resultName_, faceBlendedPtr);

    const volScalarField* nonOrthPtr =
        mesh_.lookupObjectPtr<volScalarField>(nonOrthogonalityName_);

    if (nonOrthogonality_)
    {
        if (!nonOrthPtr)
        {
            IOobject fieldHeader
            (
                nonOrthogonalityName_,
                mesh_.time().constant(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            );

            if (fieldHeader.typeHeaderOk<volScalarField>(true, true, false))
            {
                volScalarField* vfPtr(new volScalarField(fieldHeader, mesh_));
                mesh_.objectRegistry::store(vfPtr);
            }
            else
            {
                FatalErrorInFunction
                    << "Field : " << nonOrthogonalityName_ << " not found."
                    << exit(FatalError);
            }
        }
    }


    const volScalarField* faceWeightsPtr =
        mesh_.lookupObjectPtr<volScalarField>(faceWeightName_);

    if (faceWeight_)
    {
        if (!faceWeightsPtr)
        {
            IOobject fieldHeader
            (
                faceWeightName_,
                mesh_.time().constant(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            );

            if (fieldHeader.typeHeaderOk<volScalarField>(true, true, false))
            {
                volScalarField* vfPtr(new volScalarField(fieldHeader, mesh_));
                mesh_.objectRegistry::store(vfPtr);
            }
            else
            {
                FatalErrorInFunction
                    << "Field : " << faceWeightName_ << " not found."
                    << exit(FatalError);
            }
        }
    }

    const volScalarField* skewnessPtr =
        mesh_.lookupObjectPtr<volScalarField>(skewnessName_);

    if (skewness_)
    {
        if (!skewnessPtr)
        {
            IOobject fieldHeader
            (
                skewnessName_,
                mesh_.time().constant(),
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            );

            if (fieldHeader.typeHeaderOk<volScalarField>(true, true, false))
            {
                volScalarField* vfPtr(new volScalarField(fieldHeader, mesh_));
                mesh_.objectRegistry::store(vfPtr);
            }
            else
            {
                FatalErrorInFunction
                    << "Field : " << skewnessName_ << " not found."
                    << exit(FatalError);
            }
        }
    }

    if (log)
    {
        indicator_.writeOpt() = IOobject::AUTO_WRITE;
    }

    init(true);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::stabilityBlendingFactor::~stabilityBlendingFactor()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::stabilityBlendingFactor::read
(
    const dictionary& dict
)
{
    if (fieldExpression::read(dict) && writeFile::read(dict))
    {
        dict.lookup("switchNonOrtho") >> nonOrthogonality_;
        dict.lookup("switchGradCc") >> gradCc_;
        dict.lookup("switchResiduals") >> residuals_;
        dict.lookup("switchFaceWeight") >> faceWeight_;
        dict.lookup("switchSkewness") >> skewness_;
        dict.lookup("switchCo") >> Co_;

        dict.readIfPresent("maxNonOrthogonality", maxNonOrthogonality_);
        dict.readIfPresent("maxGradCc", maxGradCc_);
        dict.readIfPresent("maxResidual", maxResidual_);
        dict.readIfPresent("maxSkewness", maxSkewness_);
        dict.readIfPresent("maxFaceWeight", maxFaceWeight_);
        dict.readIfPresent("Co2", Co2_);

        dict.readIfPresent("minFaceWeight", minFaceWeight_);
        dict.readIfPresent("minNonOrthogonality", minNonOrthogonality_);
        dict.readIfPresent("minGradCc", minGradCc_);
        dict.readIfPresent("minSkewness", minSkewness_);
        dict.readIfPresent("Co1", Co1_);


        dict.readIfPresent("P", P_);
        dict.readIfPresent("I", I_);
        dict.readIfPresent("D", D_);

        tolerance_ = 0.001;
        if
        (
            dict.readIfPresent("tolerance", tolerance_)
         && (tolerance_ < 0 || tolerance_ > 1)
        )
        {
            FatalErrorInFunction
                << "tolerance must be in the range 0 to 1.  Supplied value: "
                << tolerance_ << exit(FatalError);
        }

        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::functionObjects::stabilityBlendingFactor::write()
{
    // Generate scheme statistics
    label nCellsScheme1 = 0;
    label nCellsScheme2 = 0;
    label nCellsBlended = 0;
    forAll(indicator_, celli)
    {
        scalar i = indicator_[celli];

        if (i < tolerance_)
        {
            nCellsScheme2++;
        }
        else if (i > (1 - tolerance_))
        {
            nCellsScheme1++;
        }
        else
        {
            nCellsBlended++;
        }
    }

    reduce(nCellsScheme1, sumOp<label>());
    reduce(nCellsScheme2, sumOp<label>());
    reduce(nCellsBlended, sumOp<label>());

    Log << "    scheme 1 cells :  " << nCellsScheme1 << nl
        << "    scheme 2 cells :  " << nCellsScheme2 << nl
        << "    blended cells  :  " << nCellsBlended << nl
        << endl;

    writeTime(file());

    file()
        << token::TAB << nCellsScheme1
        << token::TAB << nCellsScheme2
        << token::TAB << nCellsBlended
        << endl;

    return true;
}


// ************************************************************************* //
