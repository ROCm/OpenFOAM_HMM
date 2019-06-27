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

#include "objectiveIncompressible.H"
#include "incompressiblePrimalSolver.H"
#include "createZeroField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(objectiveIncompressible, 0);
defineRunTimeSelectionTable(objectiveIncompressible, dictionary);
addToRunTimeSelectionTable
(
    objective,
    objectiveIncompressible,
    objective
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

objectiveIncompressible::objectiveIncompressible
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& adjointSolverName,
    const word& primalSolverName
)
:
    objective(mesh, dict, adjointSolverName, primalSolverName),

    vars_
    (
        mesh.lookupObject<incompressiblePrimalSolver>(primalSolverName).
            getVars()
    ),

    // Initialize pointers to nullptr.
    // Not all of them are required for each objective function.
    // Each child should allocate whatever is needed.

    // Field adjoint Eqs
    dJdvPtr_(nullptr),
    dJdpPtr_(nullptr),
    dJdTPtr_(nullptr),
    dJdTMvar1Ptr_(nullptr),
    dJdTMvar2Ptr_(nullptr),

    // Adjoint boundary conditions
    bdJdvPtr_(nullptr),
    bdJdvnPtr_(nullptr),
    bdJdvtPtr_(nullptr),
    bdJdpPtr_(nullptr),
    bdJdTPtr_(nullptr),
    bdJdTMvar1Ptr_(nullptr),
    bdJdTMvar2Ptr_(nullptr)
{
    weight_ = dict.get<scalar>("weight");
    computeMeanFields_ = vars_.computeMeanFields();
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<objectiveIncompressible> objectiveIncompressible::New
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& adjointSolverName,
    const word& primalSolverName
)
{
    const word objectiveName = dict.dictName();
    const word modelType(dict.get<word>("type"));

    Info<< "Creating objective function : " << objectiveName
        << " of type " << modelType << endl;

    auto cstrIter = dictionaryConstructorTablePtr_->cfind(modelType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown objectiveIncompressible type " << modelType << nl << nl
            << "Valid objectiveIncompressible types are :" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<objectiveIncompressible>
    (
        cstrIter()(mesh, dict, adjointSolverName, primalSolverName)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const volVectorField& objectiveIncompressible::dJdv()
{
    if (dJdvPtr_.empty())
    {
        // If pointer is not set, set it to a zero field
        dJdvPtr_.reset
        (
            createZeroFieldPtr<vector>
            (
                mesh_,
                ("dJdv_"+type()),
                dimensionSet(0, 3, -2, 0, 0, 0, 0)
            )
        );
    }
    return dJdvPtr_();
}


const volScalarField& objectiveIncompressible::dJdp()
{
    if (dJdpPtr_.empty())
    {
        // If pointer is not set, set it to a zero field
        dJdpPtr_.reset
        (
            createZeroFieldPtr<scalar>
            (
                mesh_,
                ("dJdp_"+type()),
                dimensionSet(0, 3, -2, 0, 0, 0, 0)
            )
        );
    }
    return dJdpPtr_();
}


const volScalarField& objectiveIncompressible::dJdT()
{
    if (dJdTPtr_.empty())
    {
        // If pointer is not set, set it to a zero field
        dJdTPtr_.reset
        (
            createZeroFieldPtr<scalar>
            (
                mesh_,
                ("dJdT_"+type()),
                dimensionSet(0, 3, -2, 0, 0, 0, 0)
            )
        );
    }
    return dJdTPtr_();
}


const volScalarField& objectiveIncompressible::dJdTMvar1()
{
    if (dJdTMvar1Ptr_.empty())
    {
        // If pointer is not set, set it to a zero field
        dJdTMvar1Ptr_.reset
        (
            createZeroFieldPtr<scalar>
            (
                mesh_,
                ("dJdTMvar1_"+type()),
                dimensionSet(0, 0, -2, 0, 0, 0, 0)
            )
        );
    }
    return dJdTMvar1Ptr_();
}


const volScalarField& objectiveIncompressible::dJdTMvar2()
{
    if (dJdTMvar2Ptr_.empty())
    {
        // If pointer is not set, set it to a zero field
        dJdTMvar2Ptr_.reset
        (
            createZeroFieldPtr<scalar>
            (
                mesh_,
                ("dJdTMvar2_"+type()),
                dimensionSet(0, 3, -2, 0, 0, 0, 0)
            )
        );
    }
    return dJdTMvar2Ptr_();
}


const fvPatchVectorField& objectiveIncompressible::boundarydJdv
(
    const label patchI
)
{
    if (bdJdvPtr_.empty())
    {
        bdJdvPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
    }
    return bdJdvPtr_()[patchI];
}


const fvPatchScalarField& objectiveIncompressible::boundarydJdvn
(
    const label patchI
)
{
    if (bdJdvnPtr_.empty())
    {
        bdJdvnPtr_.reset(createZeroBoundaryPtr<scalar>(mesh_));
    }
    return bdJdvnPtr_()[patchI];
}


const fvPatchVectorField& objectiveIncompressible::boundarydJdvt
(
    const label patchI
)
{
    if (bdJdvtPtr_.empty())
    {
        bdJdvtPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
    }
    return bdJdvtPtr_()[patchI];
}


const fvPatchVectorField& objectiveIncompressible::boundarydJdp
(
    const label patchI
)
{
    if (bdJdpPtr_.empty())
    {
        bdJdpPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
    }
    return bdJdpPtr_()[patchI];
}


const fvPatchScalarField& objectiveIncompressible::boundarydJdT
(
    const label patchI
)
{
    if (bdJdTPtr_.empty())
    {
        bdJdTPtr_.reset(createZeroBoundaryPtr<scalar>(mesh_));
    }
    return bdJdTPtr_()[patchI];
}


const fvPatchScalarField& objectiveIncompressible::boundarydJdTMvar1
(
    const label patchI
)
{
    if (bdJdTMvar1Ptr_.empty())
    {
        bdJdTMvar1Ptr_.reset(createZeroBoundaryPtr<scalar>(mesh_));
    }
    return bdJdTMvar1Ptr_()[patchI];
}


const fvPatchScalarField& objectiveIncompressible::boundarydJdTMvar2
(
    const label patchI
)
{
    if (bdJdTMvar2Ptr_.empty())
    {
        bdJdTMvar2Ptr_.reset(createZeroBoundaryPtr<scalar>(mesh_));
    }
    return bdJdTMvar2Ptr_()[patchI];
}


const boundaryVectorField& objectiveIncompressible::boundarydJdv()
{
    if (bdJdvPtr_.empty())
    {
        bdJdvPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
    }
    return bdJdvPtr_();
}


const boundaryScalarField& objectiveIncompressible::boundarydJdvn()
{
    if (bdJdvnPtr_.empty())
    {
        bdJdvnPtr_.reset(createZeroBoundaryPtr<scalar>(mesh_));
    }
    return bdJdvnPtr_();
}


const boundaryVectorField& objectiveIncompressible::boundarydJdvt()
{
    if (bdJdvtPtr_.empty())
    {
        bdJdvtPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
    }
    return bdJdvtPtr_();
}


const boundaryVectorField& objectiveIncompressible::boundarydJdp()
{
    if (bdJdpPtr_.empty())
    {
        bdJdpPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
    }
    return bdJdpPtr_();
}


const boundaryScalarField& objectiveIncompressible::boundarydJdT()
{
    if (bdJdTPtr_.empty())
    {
        bdJdTPtr_.reset(createZeroBoundaryPtr<scalar>(mesh_));
    }
    return bdJdTPtr_();
}


const boundaryScalarField& objectiveIncompressible::boundarydJdTMvar1()
{
    if (bdJdTMvar1Ptr_.empty())
    {
        bdJdTMvar1Ptr_.reset(createZeroBoundaryPtr<scalar>(mesh_));
    }
    return bdJdTMvar1Ptr_();
}


const boundaryScalarField& objectiveIncompressible::boundarydJdTMvar2()
{
    if (bdJdTMvar2Ptr_.empty())
    {
        bdJdTMvar2Ptr_.reset(createZeroBoundaryPtr<scalar>(mesh_));
    }
    return bdJdTMvar2Ptr_();
}


void objectiveIncompressible::update()
{
    // Objective function value
    J();

    // Update mean values here since they might be used in the
    // subsequent functions
    update_meanValues();

    // volFields
    update_dJdv();
    update_dJdp();
    update_dJdT();
    update_dJdTMvar1();
    update_dJdTMvar2();
    update_dJdb();
    update_divDxDbMultiplier();
    update_gradDxDbMultiplier();

    // boundaryFields
    update_boundarydJdv();
    update_boundarydJdvn();
    update_boundarydJdvt();
    update_boundarydJdp();
    update_boundarydJdT();
    update_boundarydJdTMvar1();
    update_boundarydJdTMvar2();
    update_boundarydJdb();
    update_dSdbMultiplier();
    update_dndbMultiplier();
    update_dxdbMultiplier();
    update_dxdbDirectMultiplier();
    update_boundaryEdgeContribution();
    update_dJdStressMultiplier();
}


void objectiveIncompressible::write() const
{
    objective::write();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
