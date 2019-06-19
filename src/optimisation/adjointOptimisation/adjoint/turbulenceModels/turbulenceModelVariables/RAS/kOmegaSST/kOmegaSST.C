/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2007-2019 PCOpt/NTUA
                            | Copyright (C) 2013-2019 FOSS GP
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

#include "kOmegaSST.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASVariables
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(kOmegaSST, 0);
addToRunTimeSelectionTable(RASModelVariables, kOmegaSST, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kOmegaSST::kOmegaSST
(
    const fvMesh& mesh,
    const solverControl& SolverControl
)
:
    RASModelVariables(mesh, SolverControl)
{
    hasTMVar1_ = true;
    TMVar1Ptr_ = mesh_.getObjectPtr<volScalarField>("k");
    TMVar1BaseName_ = "k";

    hasTMVar2_ = true;
    TMVar2Ptr_ = mesh_.getObjectPtr<volScalarField>("omega");
    TMVar2BaseName_ = "omega";

    hasNut_ = true;
    nutPtr_ = mesh_.getObjectPtr<volScalarField>("nut");

    allocateInitValues();
    allocateMeanFields();
}


kOmegaSST::~kOmegaSST()
{
    // nullify pointers
    TMVar1Ptr_ = nullptr;
    TMVar2Ptr_ = nullptr;
    nutPtr_    = nullptr;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void kOmegaSST::correctBoundaryConditions
(
    const incompressible::turbulenceModel& turbulence
)
{
    // The presence of G is required to update the boundary value of omega
    const volVectorField& U(turbulence.U());
    const volScalarField S2(2*magSqr(symm(fvc::grad(U))));
    volScalarField G(turbulence.GName(), nutRef() * S2);
    RASModelVariables::correctBoundaryConditions(turbulence);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASVariables
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
