/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 Upstream CFD GmbH
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "SLADelta.H"
#include "wallDist.H"
#include "maxDeltaxyz.H"
#include "fvcGrad.H"
#include "fvcCurl.H"
#include "addToRunTimeSelectionTable.H"

#include "oneField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type, class WeightType = oneField>
tmp<scalarField> sumNeighbours
(
    const GeometricField<Type, fvPatchField, volMesh>& field,
    GeometricField<Type, fvPatchField, volMesh>& result,
    const WeightType& weights = oneField()
)
{
    const fvMesh& mesh = field.mesh();

    result == dimensioned<Type>(field.dimensions(), Zero);

    auto tscaling = tmp<scalarField>::New(mesh.nCells(), Zero);
    auto& scaling = tscaling.ref();

    const auto& faceNeighbour = mesh.faceNeighbour();
    const auto& faceOwner = mesh.faceOwner();

    forAll(faceNeighbour, i)
    {
        const label own = faceOwner[i];
        const label nbr = faceNeighbour[i];

        result[own] += field[nbr];
        result[nbr] += field[own];

        scaling[own] = scaling[own] + weights[nbr];
        scaling[nbr] = scaling[nbr] + weights[own];
    }

    const auto& pbm = mesh.boundaryMesh();

    for (const polyPatch& pp : pbm)
    {
        const auto& pf = field.boundaryField()[pp.index()];
        const labelUList& faceCells = pp.faceCells();

        if (pf.coupled())
        {
            const scalarField nbrField(pf.patchNeighbourField());

            forAll(faceCells, facei)
            {
                const label celli = faceCells[facei];
                result[celli] += nbrField[facei];
                scaling[celli] = scaling[celli] + weights[celli];
            }
        }
    }

    return tscaling;
}

}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace LESModels
{
    defineTypeNameAndDebug(DeltaSLA, 0);
    addToRunTimeSelectionTable(LESdelta, DeltaSLA, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::LESModels::DeltaSLA::calcDelta()
{
    const fvMesh& mesh = turbulenceModel_.mesh();

    const volVectorField& U0 = turbulenceModel_.U();
    const volVectorField vorticity(fvc::curl(U0));
    const volSymmTensorField S(symm(fvc::grad(U0)));
    const volScalarField magGradU(mag(fvc::grad(U0)));

    tmp<volScalarField> tnut = turbulenceModel_.nut();
    const auto& nut = tnut();
    tmp<volScalarField> tnu = turbulenceModel_.nu();
    const auto& nu = tnu();

    // Vortex tilting measure, VTM
    const dimensionedScalar nuMin("SMALL", nu.dimensions(), SMALL);
    const dimensionedScalar eps("SMALL", dimless/pow3(dimTime), SMALL);
    volScalarField vtm
    (
        max(scalar(1), 0.2*nu/max(nut, nuMin))
       *sqrt(6.0)*mag((S & vorticity)^vorticity)
       /max(magSqr(vorticity)*sqrt(3*magSqr(S) - sqr(tr(S))), eps)
    );
    vtm.correctBoundaryConditions();


    // Averaged VTM

    volScalarField vtmAve("vtmAve", vtm);
    tmp<scalarField> weights = sumNeighbours(vtm, vtmAve);

    // Add cell centre values
    vtmAve += vtm;

    // Weights normalisation (add 1 for cell centres)
    vtmAve.primitiveFieldRef() /= weights + 1;


    // Compute DDES shielding function
    const volScalarField& y = wallDist::New(mesh).y();
    const dimensionedScalar Ueps("eps", magGradU.dimensions(), SMALL);
    const dimensionedScalar nuEps("eps", nu.dimensions(), ROOTSMALL);
    const volScalarField rd
    (
        min
        (
            (nut + nu)/(max(magGradU, Ueps)*sqr(kappa_*y) + nuEps),
            scalar(10)
        )
    );
    const volScalarField fd(1.0 - tanh(pow(Cd1_*rd, Cd2_)));

    // Assemble delta
    const dimensionedScalar vortMin("SMALL", dimless/dimTime, SMALL);
    const volVectorField nVecVort(vorticity/(max(mag(vorticity), vortMin)));

    const cellList& cells = mesh.cells();
    const vectorField& cellCentres = mesh.cellCentres();
    const vectorField& faceCentres = mesh.faceCentres();

    forAll(cells, celli)
    {
        const labelList& cFaces = cells[celli];
        const point& cc = cellCentres[celli];
        const vector& nv = nVecVort[celli];

        scalar deltaMaxTmp = 0.0;

        for (const label facei : cFaces)
        {
            const point& fc = faceCentres[facei];
            scalar tmp = 2.0*mag(nv ^ (fc - cc));

            if (tmp > deltaMaxTmp)
            {
                deltaMaxTmp = tmp;
            }
        }

        scalar FKH =
            max
            (
                FKHmin_,
                min
                (
                    FKHmax_,
                    FKHmin_
                  + (FKHmax_ - FKHmin_)*(vtmAve[celli] - a1_)/(a2_ - a1_)
                )
            );

        if ((shielding_ == true) && (fd[celli] < (1.0 - epsilon_)))
        {
            FKH = 1;
        }

        delta_[celli] = deltaCoeff_*deltaMaxTmp*FKH;
    }

    label nD = mesh.nGeometricD();

    if (nD == 2)
    {
        WarningInFunction
            << "Case is 2D, LES is not strictly applicable" << nl
            << endl;
    }
    else if (nD != 3)
    {
        FatalErrorInFunction
            << "Case must be either 2D or 3D" << exit(FatalError);
    }

    // Handle coupled boundaries
    delta_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LESModels::DeltaSLA::DeltaSLA
(
    const word& name,
    const turbulenceModel& turbulence,
    const dictionary& dict
)
:
    LESdelta(name, turbulence),
    hmaxPtr_(nullptr),
    deltaCoeff_
    (
        dict.optionalSubDict(type() + "Coeffs").lookupOrDefault<scalar>
        (
            "deltaCoeff",
            1.035
        )
    ),
    requireUpdate_
    (
        dict.optionalSubDict(type() + "Coeffs").getOrDefault<bool>
        (
            "requireUpdate",
            true
        )
    ),
    shielding_
    (
        dict.optionalSubDict(type() + "Coeffs").getOrDefault<bool>
        (
            "shielding",
            true
        )
    ),
    FKHmin_
    (
        dict.optionalSubDict(type() + "Coeffs").getOrDefault<scalar>
        (
            "FKHmin",
            0.1
        )
    ),
    FKHmax_
    (
        dict.optionalSubDict(type() + "Coeffs").getOrDefault<scalar>
        (
            "FKHmax",
            1.0
        )
    ),
    a1_
    (
        dict.optionalSubDict(type() + "Coeffs").getOrDefault<scalar>
        (
            "a1",
            0.15
        )
    ),
    a2_
    (
        dict.optionalSubDict(type() + "Coeffs").getOrDefault<scalar>
        (
            "a2",
            0.30
        )
    ),
    epsilon_
    (
        dict.optionalSubDict(type() + "Coeffs").getOrDefault<scalar>
        (
            "epsilon",
            0.01
        )
    ),
    kappa_
    (
        dict.optionalSubDict(type() + "Coeffs").getOrDefault<scalar>
        (
            "kappa",
            0.41
        )
    ),
    Cd1_
    (
        dict.optionalSubDict(type() + "Coeffs").getOrDefault<scalar>
        (
            "Cd1",
            20
        )
    ),
    Cd2_
    (
        dict.optionalSubDict(type() + "Coeffs").getOrDefault<scalar>
        (
            "Cd2",
            3
        )
    )
{
    if (dict.optionalSubDict(type() + "Coeffs").found("hmax"))
    {
        // User-defined hmax
        hmaxPtr_ =
            LESdelta::New
            (
                IOobject::groupName("hmax", turbulence.U().group()),
                turbulence,
                dict.optionalSubDict("hmaxCoeffs"),
                "hmax"
            );
    }
    else
    {
        Info<< "Employing " << maxDeltaxyz::typeName << " for hmax" << endl;

        hmaxPtr_.reset
        (
            new maxDeltaxyz
            (
                IOobject::groupName("hmax", turbulence.U().group()),
                turbulence,
                dict.optionalSubDict("hmaxCoeffs")
            )
        );
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::LESModels::DeltaSLA::read(const dictionary& dict)
{
    const dictionary& coeffsDict(dict.optionalSubDict(type() + "Coeffs"));

    coeffsDict.readIfPresent<scalar>("deltaCoeff", deltaCoeff_);
    coeffsDict.readIfPresent<bool>("requireUpdate", requireUpdate_);
    coeffsDict.readIfPresent<scalar>("FKHmin", FKHmin_);
    coeffsDict.readIfPresent<scalar>("FKHmax", FKHmax_);
    coeffsDict.readIfPresent<scalar>("a1", a1_);
    coeffsDict.readIfPresent<scalar>("a2", a2_);
    coeffsDict.readIfPresent<scalar>("epsilon", epsilon_);
    coeffsDict.readIfPresent<scalar>("kappa", kappa_);
    coeffsDict.readIfPresent<scalar>("Cd1", Cd1_);
    coeffsDict.readIfPresent<scalar>("Cd2", Cd2_);

    calcDelta();
}


void Foam::LESModels::DeltaSLA::correct()
{
    if (turbulenceModel_.mesh().changing() && requireUpdate_)
    {
        hmaxPtr_->correct();
    }

    calcDelta();
}


// ************************************************************************* //
