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

#include "adjointBoundaryCondition.H"
#include "ATCUaGradU.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(adjointBoundaryCondition, 0);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool adjointBoundaryCondition::addATCUaGradUTerm()
{
    if (addATCUaGradUTerm_.empty())
    {
        addATCUaGradUTerm_.reset(new bool(isA<ATCUaGradU>(getATC())));
    }
    return addATCUaGradUTerm_();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

adjointBoundaryCondition::adjointBoundaryCondition
(
    const adjointBoundaryCondition& adjointBC
)
:
    patch_(adjointBC.patch_),
    managerName_(adjointBC.managerName_),
    adjointSolverName_(adjointBC.adjointSolverName_),
    simulationType_(adjointBC.simulationType_),
    boundaryContrPtr_
    (
        boundaryAdjointContribution::New
        (
            adjointBC.managerName_,
            adjointBC.adjointSolverName_,
            adjointBC.simulationType_,
            adjointBC.patch_
        )
    ),
    addATCUaGradUTerm_(adjointBC.addATCUaGradUTerm_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const word& adjointBoundaryCondition::objectiveManagerName() const
{
    return managerName_;
}


const word& adjointBoundaryCondition::adjointSolverName() const
{
    return adjointSolverName_;
}


const word& adjointBoundaryCondition::simulationType() const
{
    return simulationType_;
}


void adjointBoundaryCondition::setBoundaryContributionPtr()
{
    // Note:
    // Check whether there is an objectiveFunctionManager object in the registry
    // Necessary for decomposePar if the libadjoint is loaded
    // through controlDict. A nicer way should be found
    const fvMesh& meshRef = patch_.boundaryMesh().mesh();
    if (meshRef.foundObject<regIOobject>(managerName_))
    {
        boundaryContrPtr_.reset
        (
            boundaryAdjointContribution::New
            (
                managerName_,
                adjointSolverName_,
                simulationType_,
                patch_
            ).ptr()
        );
    }
    else
    {
        WarningInFunction
            << "No objectiveManager " << managerName_ << " available." << nl
            << "Setting boundaryAdjointContributionPtr to nullptr. " << nl
            << "OK for decomposePar."
            << endl;
    }
}


boundaryAdjointContribution&
adjointBoundaryCondition::getBoundaryAdjContribution()
{
    return boundaryContrPtr_();
}


const ATCModel& adjointBoundaryCondition::getATC() const
{
    return
        patch_.boundaryMesh().mesh().lookupObject<ATCModel>
        (
            "ATCModel" + adjointSolverName_
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
