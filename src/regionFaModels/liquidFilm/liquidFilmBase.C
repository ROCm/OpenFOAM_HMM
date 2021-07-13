/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "liquidFilmBase.H"
#include "faMesh.H"
#include "faCFD.H"
#include "uniformDimensionedFields.H"
#include "gravityMeshObject.H"
#include "movingWallVelocityFvPatchVectorField.H"
#include "turbulentFluidThermoModel.H"
#include "turbulentTransportModel.H"
#include "calculatedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(liquidFilmBase, 0);

defineRunTimeSelectionTable(liquidFilmBase, dictionary);

const Foam::word liquidFilmName("liquidFilm");

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //


bool liquidFilmBase::read(const dictionary& dict)
{
    regionFaModel::read(dict);
    if (active_)
    {
        const dictionary& solution = this->solution().subDict("PISO");
        solution.readEntry("momentumPredictor", momentumPredictor_);
        solution.readIfPresent("nOuterCorr", nOuterCorr_);
        solution.readEntry("nCorr", nCorr_);
        solution.readEntry("nNonOrthCorr", nNonOrthCorr_);
    }
    return true;
}


scalar liquidFilmBase::CourantNumber() const
{
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;
    scalar velMag = 0.0;

    if (regionMesh().nInternalEdges())
    {
        edgeScalarField SfUfbyDelta
        (
            regionMesh().edgeInterpolation::deltaCoeffs()*mag(phif_)
        );

        CoNum = max(SfUfbyDelta/regionMesh().magLe())
            .value()*time().deltaT().value();

        meanCoNum = (sum(SfUfbyDelta)/sum(regionMesh().magLe()))
            .value()*time().deltaT().value();

        velMag = max(mag(phif_)/regionMesh().magLe()).value();
    }

    Info<< "Film Courant Number mean: " << meanCoNum
        << " max: " << CoNum
        << " Film velocity magnitude: " << velMag << endl;

    return CoNum;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


liquidFilmBase::liquidFilmBase
(
    const word& modelType,
    const fvPatch& p,
    const dictionary& dict
)
:
    regionFaModel(p, liquidFilmName, modelType, dict, true),

    momentumPredictor_
    (
        this->solution().subDict("PISO").lookup("momentumPredictor")
    ),
    nOuterCorr_
    (
        this->solution().subDict("PISO").lookupOrDefault("nOuterCorr", 1)
    ),
    nCorr_(this->solution().subDict("PISO").get<label>("nCorr")),
    nNonOrthCorr_
    (
        this->solution().subDict("PISO").get<label>("nNonOrthCorr")
    ),

    h0_("h0", dimLength, SMALL),

    UName_(dict.get<word>("U")),

    pName_(dict.lookupOrDefault<word>("p",  word::null)),

    pRef_(dict.get<scalar>("pRef")),

    h_
    (
        IOobject
        (
            "hf_" + regionName_,
            primaryMesh().time().timeName(),
            primaryMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh()
    ),

    Uf_
    (
        IOobject
        (
            "Uf_" + regionName_,
            primaryMesh().time().timeName(),
            primaryMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh()
    ),
    pf_
    (
        IOobject
        (
            "pf_" + regionName_,
            primaryMesh().time().timeName(),
            primaryMesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimPressure, Zero)
    ),
    ppf_
    (
        IOobject
        (
            "ppf_" + regionName_,
            primaryMesh().time().timeName(),
            primaryMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimPressure, Zero)
    ),
    phif_
    (
        IOobject
        (
            "phif_" + regionName_,
            primaryMesh().time().timeName(),
            primaryMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fac::interpolate(Uf_) & regionMesh().Le()
    ),

    phi2s_
    (
        IOobject
        (
            "phi2s_" + regionName_,
            primaryMesh().time().timeName(),
            primaryMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fac::interpolate(h_*Uf_) & regionMesh().Le()
    ),

    Tf_
    (
        IOobject
        (
            "Tf_" + regionName_,
            primaryMesh().time().timeName(),
            primaryMesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        regionMesh()
    ),

    gn_
    (
        IOobject
        (
            "gn",
            primaryMesh().time().timeName(),
            primaryMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        regionMesh(),
        dimensionedScalar(dimAcceleration, Zero)
    ),

    massSource_
    (
        IOobject
        (
            "massSource",
            primaryMesh().time().timeName(),
            primaryMesh()
        ),
        primaryMesh(),
        dimensionedScalar(dimMass, Zero),
        calculatedFvPatchField<scalar>::typeName
    ),

    momentumSource_
    (
        IOobject
        (
            "momentumSource",
            primaryMesh().time().timeName(),
            primaryMesh()
        ),
        primaryMesh(),
        dimensionedVector(dimPressure, Zero),
        calculatedFvPatchField<vector>::typeName
    ),

    pnSource_
    (
        IOobject
        (
            "pnSource",
            primaryMesh().time().timeName(),
            primaryMesh()
        ),
        primaryMesh(),
        dimensionedScalar(dimPressure, Zero),
        calculatedFvPatchField<scalar>::typeName
    ),

    energySource_
    (
        IOobject
        (
            "energySource",
            primaryMesh().time().timeName(),
            primaryMesh()
        ),
        primaryMesh(),
        dimensionedScalar(dimEnergy, Zero),
        calculatedFvPatchField<scalar>::typeName
    ),

    addedMassTotal_(0),

    faOptions_(Foam::fa::options::New(p))
{
    const uniformDimensionedVectorField& g =meshObjects::gravity::New(time());

    gn_ = (g & regionMesh().faceAreaNormals());


    if (!faOptions_.optionList::size())
    {
        Info << "No finite area options present" << endl;
    }
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

liquidFilmBase::~liquidFilmBase()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::areaVectorField> liquidFilmBase::Uw() const
{
    tmp<areaVectorField> tUw
    (
        new areaVectorField
        (
            IOobject
            (
                "tUw",
                primaryMesh().time().timeName(),
                primaryMesh()
            ),
            regionMesh(),
            dimensionedVector(dimVelocity, Zero)
        )
    );

    areaVectorField& Uw = tUw.ref();

    const polyPatch& pp = primaryMesh().boundaryMesh()[patch_.index()];

    if
    (
        primaryMesh().moving()
     && isA<movingWallVelocityFvPatchVectorField>(pp)
    )
    {
        const movingWallVelocityFvPatchVectorField& wpp =
            refCast<const movingWallVelocityFvPatchVectorField>(pp);

        tmp<vectorField> tUwall = wpp.Uwall();

         // Map Up to area
        tmp<vectorField> UsWall = vsmPtr_->mapToSurface(tUwall());

        const vectorField& nHat =
            regionMesh().faceAreaNormals().internalField();

        Uw.primitiveFieldRef() = UsWall() - nHat*(UsWall() & nHat);
    }

    return tUw;
}


Foam::tmp<Foam::areaVectorField> liquidFilmBase::Us() const
{
    tmp<areaVectorField> tUs
    (
        new areaVectorField
        (
            IOobject
            (
                "tUs",
                primaryMesh().time().timeName(),
                primaryMesh()
            ),
            regionMesh(),
            dimensionedVector(dimVelocity, Zero)
        )
    );

    // Uf quadratic profile
    tUs.ref() = Foam::sqrt(2.0)*Uf_;

    return tUs;
}


Foam::tmp<Foam::areaVectorField> liquidFilmBase::Up() const
{
    const label patchi = patch_.index();

    const volVectorField& Uprimary =
        primaryMesh().lookupObject<volVectorField>(UName_);

    const fvPatchVectorField& Uw = Uprimary.boundaryField()[patchi];

    tmp<areaVectorField> tUp
    (
        new areaVectorField
        (
            IOobject
            (
                "tUp",
                primaryMesh().time().timeName(),
                primaryMesh()
            ),
            regionMesh(),
            dimensionedVector(dimVelocity, Zero)
        )
    );

    areaVectorField& Up = tUp.ref();
/*
    typedef compressible::turbulenceModel cmpTurbModelType;
    typedef incompressible::turbulenceModel incmpTurbModelType;

    word turbName
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            Up_.internalField().group()
        )
    );

    scalarField nu(patch_.size(), Zero);
    if (primaryMesh().foundObject<cmpTurbModelType>(turbName))
    {
        const cmpTurbModelType& turbModel =
            primaryMesh().lookupObject<cmpTurbModelType>(turbName);

        //const basicThermo& thermo = turbModel.transport();

        nu = turbModel.nu(patchi);
    }
    if (primaryMesh().foundObject<incmpTurbModelType>(turbName))
    {
        const incmpTurbModelType& turbModel =
            primaryMesh().lookupObject<incmpTurbModelType>(turbName);

        nu = turbModel.nu(patchi);
    }
*/
    scalarField hp(patch_.size(), Zero);

    // map areas h to hp on primary
    vsmPtr_->mapToField(h_, hp);

    const vectorField& nHat = regionMesh().faceAreaNormals().internalField();

    // U tangential on primary
    const vectorField Ust(-Uw.snGrad()*hp);

    // Map U tang to surface
    Up.primitiveFieldRef() = vsmPtr_->mapToSurface(Ust);

    // U tangent on surface
    Up.primitiveFieldRef() -= nHat*(Up.primitiveField() & nHat);

    return tUp;
}

tmp<areaScalarField> liquidFilmBase::pg() const
{
    tmp<areaScalarField> tpg
    (
        new areaScalarField
        (
            IOobject
            (
                "tpg",
                primaryMesh().time().timeName(),
                primaryMesh()
            ),
            regionMesh(),
            dimensionedScalar(dimPressure, Zero)
        )
    );

    areaScalarField& pfg = tpg.ref();

    if (pName_ != word::null)
    {
        const volScalarField& pp =
            primaryMesh().lookupObject<volScalarField>(pName_);

        const volScalarField::Boundary& pw = pp.boundaryField();

        pfg.primitiveFieldRef() = vsmPtr_->mapInternalToSurface<scalar>(pw)();
    }

    return tpg;
}


void liquidFilmBase::addSources
(
    const label patchi,
    const label facei,
    const scalar massSource,
    const vector& momentumSource,
    const scalar pressureSource,
    const scalar energySource
)
{

    massSource_.boundaryFieldRef()[patchi][facei] += massSource;

    pnSource_.boundaryFieldRef()[patchi][facei] += pressureSource;

    momentumSource_.boundaryFieldRef()[patchi][facei] += momentumSource;

    addedMassTotal_ += massSource;

    if (debug)
    {
        InfoInFunction
            << "\nSurface film: " << type() << ": adding to film source:" << nl
            << "    mass     = " << massSource << nl
            << "    momentum = " << momentumSource << nl
            << "    pressure = " << pressureSource << endl;
    }
}


Foam::fa::options& liquidFilmBase::faOptions()
{
     return faOptions_;
}


const areaVectorField& liquidFilmBase::Uf() const
{
     return Uf_;
}


const areaScalarField& liquidFilmBase::Tf() const
{
     return Tf_;
}

const areaScalarField& liquidFilmBase::gn() const
{
     return gn_;
}

const areaScalarField& liquidFilmBase::h() const
{
     return h_;
}

const dimensionedScalar& liquidFilmBase::h0() const
{
     return h0_;
}

const regionFaModel& liquidFilmBase::region() const
{
    return *this;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
