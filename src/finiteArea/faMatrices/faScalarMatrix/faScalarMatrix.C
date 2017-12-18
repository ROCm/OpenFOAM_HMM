/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2016-2017 Wikki Ltd
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

Description
     Finite-Area scalar matrix member functions and operators

\*---------------------------------------------------------------------------*/

#include "faScalarMatrix.H"
#include "zeroGradientFaPatchFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
void Foam::faMatrix<Foam::scalar>::setComponentReference
(
    const label patchI,
    const label edgeI,
    const direction,
    const scalar value
)
{
    const labelUList& faceLabels = psi_.mesh().boundary()[patchI].edgeFaces();

    internalCoeffs_[patchI][edgeI] += diag()[faceLabels[edgeI]];

    boundaryCoeffs_[patchI][edgeI] = value;
}


template<>
Foam::solverPerformance Foam::faMatrix<Foam::scalar>::solve
(
    const dictionary& solverControls
)
{
    DebugInFunction
        << "solving faMatrix<scalar>"
        << endl;

    GeometricField<scalar, faPatchField, areaMesh>& psi =
        const_cast<GeometricField<scalar, faPatchField, areaMesh>&>(psi_);

    scalarField saveDiag(diag());
    addBoundaryDiag(diag(), 0);

    scalarField totalSource(source_);
    addBoundarySource(totalSource, 0);

    // Solver call
    solverPerformance solverPerf = lduMatrix::solver::New
    (
        psi_.name(),
        *this,
        boundaryCoeffs_,
        internalCoeffs_,
        psi_.boundaryField().scalarInterfaces(),
        solverControls
    )->solve(psi.ref(), totalSource);

    if (solverPerformance::debug)
    {
        solverPerf.print(Info);
    }

    diag() = saveDiag;

    psi.correctBoundaryConditions();

    psi.mesh().setSolverPerformance(psi.name(), solverPerf);

    return solverPerf;
}


template<>
Foam::tmp<Foam::scalarField> Foam::faMatrix<Foam::scalar>::residual() const
{
    scalarField boundaryDiag(psi_.size(), 0.0);
    addBoundaryDiag(boundaryDiag, 0);

    tmp<scalarField> tres
    (
        lduMatrix::residual
        (
            psi_.internalField(),
            source_ - boundaryDiag*psi_.internalField(),
            boundaryCoeffs_,
            psi_.boundaryField().scalarInterfaces(),
            0
        )
    );

    addBoundarySource(tres.ref());

    return tres;
}


template<>
Foam::tmp<Foam::areaScalarField> Foam::faMatrix<Foam::scalar>::H() const
{
    tmp<areaScalarField> tHphi
    (
        new areaScalarField
        (
            IOobject
            (
                "H("+psi_.name()+')',
                psi_.instance(),
                psi_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            psi_.mesh(),
            dimensions_/dimArea,
            zeroGradientFaPatchScalarField::typeName
        )
    );
    areaScalarField& Hphi = tHphi.ref();

    Hphi.primitiveFieldRef() = (lduMatrix::H(psi_.primitiveField()) + source_);
    addBoundarySource(Hphi.primitiveFieldRef());

    Hphi.ref() /= psi_.mesh().S();
    Hphi.correctBoundaryConditions();

    return tHphi;
}


// ************************************************************************* //
