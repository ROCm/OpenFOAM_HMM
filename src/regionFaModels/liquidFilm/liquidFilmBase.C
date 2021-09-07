/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2021 OpenCFD Ltd.
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
        this->solution().subDict("PIMPLE").get<bool>("momentumPredictor")
    ),
    nOuterCorr_
    (
        this->solution().subDict("PIMPLE").get<label>("nOuterCorr")
    ),
    nCorr_(this->solution().subDict("PIMPLE").get<label>("nCorr")),
    nFilmCorr_
    (
        this->solution().subDict("PIMPLE").get<label>("nFilmCorr")
    ),

    h0_("h0", dimLength, 1e-7, dict),

    deltaWet_("deltaWet", dimLength, 1e-4, dict),

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

    g_(meshObjects::gravity::New(primaryMesh().time())),

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
    const areaVectorField& ns = regionMesh().faceAreaNormals();

    gn_ = g_ & ns;

    if (!faOptions_.optionList::size())
    {
        Info << "No finite area options present" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

liquidFilmBase::~liquidFilmBase()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

scalar liquidFilmBase::CourantNumber() const
{
    scalar CoNum = 0.0;
    scalar velMag = 0.0;

    edgeScalarField SfUfbyDelta
    (
        regionMesh().edgeInterpolation::deltaCoeffs()*mag(phif_)
    );

    CoNum =
        max(SfUfbyDelta/regionMesh().magLe()).value()*time().deltaT().value();

    velMag = max(mag(phif_)/regionMesh().magLe()).value();

    reduce(CoNum, maxOp<scalar>());
    reduce(velMag, maxOp<scalar>());

    Info<< "Max film Courant Number: " << CoNum
        << " Film velocity magnitude: " << velMag << endl;

    return CoNum;
}


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

        volScalarField::Boundary& pw =
            const_cast<volScalarField::Boundary&>(pp.boundaryField());

        //pw -= pRef_;

        pfg.primitiveFieldRef() = vsmPtr_->mapInternalToSurface<scalar>(pw)();
    }
    return tpg;
}


tmp<areaScalarField> liquidFilmBase::alpha() const
{
    tmp<areaScalarField> talpha
    (
        new areaScalarField
        (
            IOobject
            (
                "talpha",
                primaryMesh().time().timeName(),
                primaryMesh()
            ),
            regionMesh(),
            dimensionedScalar(dimless, Zero)
        )
    );
    areaScalarField& alpha = talpha.ref();

    alpha = pos0(h_ - deltaWet_);

    return talpha;
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
}


void liquidFilmBase::preEvolveRegion()
{
    regionFaModel::preEvolveRegion();
}


void liquidFilmBase::postEvolveRegion()
{
    if (debug && primaryMesh().time().writeTime())
    {
        massSource_.write();
        pnSource_.write();
        momentumSource_.write();
    }

    massSource_.boundaryFieldRef() = Zero;
    pnSource_.boundaryFieldRef() = Zero;
    momentumSource_.boundaryFieldRef() = Zero;

    regionFaModel::postEvolveRegion();
}


Foam::fa::options& liquidFilmBase::faOptions()
{
     return faOptions_;
}


const areaVectorField& liquidFilmBase::Uf() const
{
     return Uf_;
}


const areaScalarField& liquidFilmBase::gn() const
{
     return gn_;
}


const uniformDimensionedVectorField& liquidFilmBase::g() const
{
    return g_;
}


const areaScalarField& liquidFilmBase::h() const
{
     return h_;
}


const edgeScalarField& liquidFilmBase::phif() const
{
    return phif_;
}


const edgeScalarField& liquidFilmBase::phi2s() const
{
    return phi2s_;
}


const dimensionedScalar& liquidFilmBase::h0() const
{
     return h0_;
}


const regionFaModel& liquidFilmBase::region() const
{
    return *this;
}


scalar liquidFilmBase::pRef()
{
    return pRef_;
}


word liquidFilmBase::UName() const
{
    return UName_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
