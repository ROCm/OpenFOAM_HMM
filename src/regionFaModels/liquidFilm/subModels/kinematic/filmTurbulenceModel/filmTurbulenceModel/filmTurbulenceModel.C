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

#include "filmTurbulenceModel.H"
#include "gravityMeshObject.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(filmTurbulenceModel, 0);
defineRunTimeSelectionTable(filmTurbulenceModel, dictionary);

const Enum
<
    filmTurbulenceModel::frictionMethodType
>
filmTurbulenceModel::frictionMethodTypeNames_
{
    { frictionMethodType::mquadraticProfile, "quadraticProfile" },
    { frictionMethodType::mlinearProfile, "linearProfile" },
    { frictionMethodType::mDarcyWeisbach, "DarcyWeisbach" },
    { frictionMethodType::mManningStrickler, "ManningStrickler" }
};


const Enum
<
    filmTurbulenceModel::shearMethodType
>
filmTurbulenceModel::shearMethodTypeNames_
{
    { shearMethodType::msimple, "simple" },
    { shearMethodType::mwallFunction, "wallFunction" }
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

filmTurbulenceModel::filmTurbulenceModel
(
    const word& modelType,
    liquidFilmBase& film,
    const dictionary& dict
)
:
    film_(film),
    dict_(dict.subDict(modelType + "Coeffs")),
    method_(frictionMethodTypeNames_.get("friction", dict_)),
    shearMethod_(shearMethodTypeNames_.get("shearStress", dict_)),
    rhoName_(dict_.getOrDefault<word>("rho", "rho")),
    rhoRef_(VGREAT)
{
    if (rhoName_ == "rhoInf")
    {
        rhoRef_ = dict_.get<scalar>("rhoInf");
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const liquidFilmBase& filmTurbulenceModel::film() const
{
    return film_;
}


tmp<areaScalarField> filmTurbulenceModel::Cw() const
{
    auto tCw =
        tmp<areaScalarField>::New
        (
            IOobject
            (
                "tCw",
                film_.primaryMesh().time().timeName(),
                film_.primaryMesh()
            ),
            film_.regionMesh(),
            dimensionedScalar(dimVelocity)
        );
    auto& Cw = tCw.ref();


    switch (method_)
    {
        case mquadraticProfile:
        {
            const scalarField& mu = film_.mu().primitiveField();
            const scalarField& h = film_.h().primitiveField();
            const scalarField& rho = film_.rho().primitiveField();

            const scalar h0 = film_.h0().value();

            Cw.primitiveFieldRef() = 3*mu/((h + h0)*rho);
            Cw.min(5000.0);

            break;
        }
        case mlinearProfile:
        {
            const scalarField& mu = film_.mu().primitiveField();
            const scalarField& h = film_.h().primitiveField();
            const scalarField& rho = film_.rho().primitiveField();

            const scalar h0 = film_.h0().value();

            Cw.primitiveFieldRef() = 2*mu/((h + h0)*rho);

            break;
        }
        case mDarcyWeisbach:
        {
            const uniformDimensionedVectorField& g =
                meshObjects::gravity::New(film_.primaryMesh().time());
            const vectorField& Uf = film_.Uf().primitiveField();
            const scalarField& rho = film_.rho().primitiveField();

            const scalar Cf = dict_.get<scalar>("DarcyWeisbach");

            Cw.primitiveFieldRef() = Cf*mag(g.value())*mag(Uf)/rho;

            break;
        }
        case mManningStrickler:
        {
            const uniformDimensionedVectorField& g =
                meshObjects::gravity::New(film_.primaryMesh().time());

            const vectorField& Uf = film_.Uf().primitiveField();
            const scalarField& h = film_.h().primitiveField();
            const scalar h0 = film_.h0().value();

            const scalar n = dict_.get<scalar>("n");

            Cw.primitiveFieldRef() =
                sqr(n)*mag(g.value())*mag(Uf)/cbrt(h + h0);

            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unimplemented method "
                << frictionMethodTypeNames_[method_] << nl
                << "Please set 'frictionMethod' to one of "
                << flatOutput(frictionMethodTypeNames_.sortedToc()) << nl
                << exit(FatalError);

            break;
        }
    }

    return tCw;
}


tmp<faVectorMatrix> filmTurbulenceModel::primaryRegionFriction
(
    areaVectorField& U
) const
{
    tmp<faVectorMatrix> tshearStress
    (
        new faVectorMatrix(U, sqr(U.dimensions())*sqr(dimLength))
    );

    switch (shearMethod_)
    {
        case msimple:
        {
            tmp<areaVectorField> Up = film_.Up();

            const dimensionedScalar Cf
            (
                "Cf",
                dimVelocity,
                dict_.get<scalar>("Cf")
            );

            tshearStress.ref() += - fam::Sp(Cf, U) + Cf*Up();

            break;
        }
        case mwallFunction:
        {
            tmp<volSymmTensorField> tdevRhoReff = devRhoReff();

            const volSymmTensorField::Boundary& devRhoReffb
                = tdevRhoReff().boundaryField();

            const label patchi = film_.patchID();

            const surfaceVectorField::Boundary& Sfb =
                film_.primaryMesh().Sf().boundaryField();

            vectorField fT(Sfb[patchi] & devRhoReffb[patchi]);

            const vectorField& nHat =
                film_.regionMesh().faceAreaNormals().internalField();

            // Substract normal component
            fT -= nHat*(fT & nHat);

            auto taForce = tmp<areaVectorField>::New
            (
                IOobject
                (
                    "taForce",
                    film_.primaryMesh().time().timeName(),
                    film_.primaryMesh()
                ),
                film_.regionMesh(),
                dimensionedVector(sqr(dimVelocity), Zero)
            );
            vectorField& aForce = taForce.ref().primitiveFieldRef();

            // Map ft to surface
            const vectorField afT(film_.vsm().mapToSurface(fT));

            const DimensionedField<scalar, areaMesh>& magSf =
                film_.regionMesh().S();

            aForce = afT/(film_.rho().primitiveField()*magSf);

            tshearStress.ref() += taForce();

            if (film_.regionMesh().time().writeTime())
            {
                taForce().write();
            }

            break;
        }
    }

    return tshearStress;
}


tmp<Foam::volSymmTensorField> filmTurbulenceModel::devRhoReff() const
{
    typedef compressible::turbulenceModel cmpTurbModel;
    typedef incompressible::turbulenceModel icoTurbModel;

    const fvMesh& m = film_.primaryMesh();

    const auto& U = m.lookupObject<volVectorField>(film_.UName());

    if (m.foundObject<cmpTurbModel>(cmpTurbModel::propertiesName))
    {
        const auto& turb =
            m.lookupObject<cmpTurbModel>(cmpTurbModel::propertiesName);

        return turb.devRhoReff();
    }
    else if (m.foundObject<icoTurbModel>(icoTurbModel::propertiesName))
    {
        const auto& turb =
            m.lookupObject<icoTurbModel>(icoTurbModel::propertiesName);

        return rho()*turb.devReff();
    }
    else if (m.foundObject<fluidThermo>(fluidThermo::dictName))
    {
        const auto& thermo =
            m.lookupObject<fluidThermo>(fluidThermo::dictName);

        return -thermo.mu()*dev(twoSymm(fvc::grad(U)));
    }
    else if (m.foundObject<transportModel>("transportProperties"))
    {
        const auto& laminarT =
            m.lookupObject<transportModel>("transportProperties");

        return -rho()*laminarT.nu()*dev(twoSymm(fvc::grad(U)));
    }
    else if (m.foundObject<dictionary>("transportProperties"))
    {
        const auto& transportProperties =
            m.lookupObject<dictionary>("transportProperties");

        const dimensionedScalar nu("nu", dimViscosity, transportProperties);

        return -rho()*nu*dev(twoSymm(fvc::grad(U)));
    }
    else
    {
        FatalErrorInFunction
            << "No valid model for viscous stress calculation"
            << exit(FatalError);

        return volSymmTensorField::null();
    }
}


tmp<Foam::volScalarField> filmTurbulenceModel::rho() const
{
    const fvMesh& m = film_.primaryMesh();
    if (rhoName_ == "rhoInf")
    {
        return tmp<volScalarField>::New
        (
            IOobject
            (
                "rho",
                m.time().timeName(),
                m
            ),
            m,
            dimensionedScalar(dimDensity, rhoRef_)
        );
    }

    return m.lookupObject<volScalarField>(rhoName_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
