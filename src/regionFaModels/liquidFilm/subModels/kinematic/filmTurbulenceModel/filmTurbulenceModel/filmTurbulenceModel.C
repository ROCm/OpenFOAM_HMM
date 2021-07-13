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

#include "filmTurbulenceModel.H"
#include "gravityMeshObject.H"

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
    method_(frictionMethodTypeNames_.get("friction", dict_))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

filmTurbulenceModel::~filmTurbulenceModel()
{}


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
            const scalar h0 = film_.h0().value();
            Cw.primitiveFieldRef() = 3*mu/(h + h0);
            Cw.min(5000.0);
            break;
        }
        case mlinearProfile:
        {
            const scalarField& mu = film_.mu().primitiveField();
            const scalarField& h = film_.h().primitiveField();
            const scalar h0 = film_.h0().value();
            Cw.primitiveFieldRef() = 2*mu/(h + h0);
            break;
        }
        case mDarcyWeisbach:
        {
            const uniformDimensionedVectorField& g =
                meshObjects::gravity::New(film_.primaryMesh().time());

            const vectorField& Uf = film_.Uf().primitiveField();

            scalar Cf = dict_.get<scalar>("Cf");
            Cw.primitiveFieldRef() = Cf*mag(g.value())*mag(Uf);
            break;
        }
        case mManningStrickler:
        {
            const uniformDimensionedVectorField& g =
                meshObjects::gravity::New(film_.primaryMesh().time());

            const scalar n = dict_.get<scalar>("n");

            const vectorField& Uf = film_.Uf().primitiveField();
            const scalarField& h = film_.h().primitiveField();
            const scalar h0 = film_.h0().value();

            Cw.primitiveFieldRef() =
                sqr(n)*mag(g.value())*mag(Uf)/cbrt(h+h0);
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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
