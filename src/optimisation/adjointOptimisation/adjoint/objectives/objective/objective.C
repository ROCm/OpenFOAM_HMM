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

#include "objective.H"
#include "createZeroField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(objective, 0);
defineRunTimeSelectionTable(objective, objective);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

void objective::makeFolder()
{
    if (Pstream::master())
    {
        const Time& time = mesh_.time();
        objFunctionFolder_ =
            time.globalPath()/"optimisation"/type()/time.timeName();

        mkDir(objFunctionFolder_);
    }
}


void objective::setObjectiveFilePtr() const
{
    objFunctionFilePtr_.reset
    (
        new OFstream(objFunctionFolder_/objectiveName_ + adjointSolverName_)
    );
}


void objective::setInstantValueFilePtr() const
{
    instantValueFilePtr_.reset
    (
        new OFstream
        (
            objFunctionFolder_/objectiveName_ + "Instant" + adjointSolverName_
        )
    );
}


void objective::setMeanValueFilePtr() const
{
    meanValueFilePtr_.reset
    (
        new OFstream
        (
            objFunctionFolder_/objectiveName_ + "Mean" + adjointSolverName_
        )
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

const dictionary& objective::dict() const
{
    return dict_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

objective::objective
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& adjointSolverName,
    const word& primalSolverName
)
:
    mesh_(mesh),
    dict_(dict),
    adjointSolverName_(adjointSolverName),
    primalSolverName_(primalSolverName),
    objectiveName_(dict.dictName()),
    computeMeanFields_(false), // is reset in derived classes

    J_(Zero),
    JMean_(Zero),
    weight_(Zero),

    // Initialize pointers to nullptr.
    // Not all of them are required for each objective function.
    // Each child should allocate whatever is needed.

    dJdbPtr_(nullptr),
    bdJdbPtr_(nullptr),
    bdSdbMultPtr_(nullptr),
    bdndbMultPtr_(nullptr),
    bdxdbMultPtr_(nullptr),
    bdxdbDirectMultPtr_(nullptr),
    bEdgeContribution_(nullptr),
    bdJdStressPtr_(nullptr),
    divDxDbMultPtr_(nullptr),
    gradDxDbMultPtr_(nullptr),

    objFunctionFolder_("word"),
    objFunctionFilePtr_(nullptr),
    instantValueFilePtr_(nullptr),
    meanValueFilePtr_(nullptr)
{
    makeFolder();
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

autoPtr<objective> objective::New
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& objectiveType,
    const word& adjointSolverName,
    const word& primalSolverName
)
{
    auto cstrIter = objectiveConstructorTablePtr_->cfind(objectiveType);

    if (!cstrIter.found())
    {
        FatalErrorInFunction
            << "Unknown objective type " << objectiveType << nl << nl
            << "Valid types are :" << nl
            << objectiveConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<objective>
    (
        cstrIter()
        (
            mesh,
            dict,
            adjointSolverName,
            primalSolverName
        )
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool objective::readDict(const dictionary& dict)
{
    dict_ = dict;
    return true;
}


void objective::updateNormalizationFactor()
{
    // Does nothing in base
}


void objective::accumulateJMean(solverControl& solverControl)
{
    if (solverControl.doAverageIter())
    {
        const label iAverageIter = solverControl.averageIter();
        if (iAverageIter == 0)
        {
            JMean_ = Zero;
        }
        scalar avIter(iAverageIter);
        scalar oneOverItP1 = 1./(avIter + 1);
        scalar mult = avIter*oneOverItP1;
        JMean_ = JMean_*mult + J_*oneOverItP1;
    }
}


scalar objective::weight() const
{
    return weight_;
}


const volScalarField& objective::dJdb()
{
    if (dJdbPtr_.empty())
    {
        // If pointer is not set, set it to a zero field
        dJdbPtr_.reset
        (
            createZeroFieldPtr<scalar>
            (
                mesh_,
                ("dJdb_" + objectiveName_),
                dimensionSet(0, 5, -2, 0, 0, 0, 0)
            )
        );
    }

    return dJdbPtr_();
}


const fvPatchVectorField& objective::boundarydJdb(const label patchI)
{
    if (bdJdbPtr_.empty())
    {
        bdJdbPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
    }
    return bdJdbPtr_()[patchI];
}


const fvPatchVectorField& objective::dSdbMultiplier(const label patchI)
{
    if (bdSdbMultPtr_.empty())
    {
        bdSdbMultPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
    }
    return bdSdbMultPtr_()[patchI];
}


const fvPatchVectorField& objective::dndbMultiplier(const label patchI)
{
    if (bdndbMultPtr_.empty())
    {
        bdndbMultPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
    }
    return bdndbMultPtr_()[patchI];
}


const fvPatchVectorField& objective::dxdbMultiplier(const label patchI)
{
    if (bdxdbMultPtr_.empty())
    {
        bdxdbMultPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
    }
    return bdxdbMultPtr_()[patchI];
}


const fvPatchVectorField& objective::dxdbDirectMultiplier(const label patchI)
{
    if (bdxdbDirectMultPtr_.empty())
    {
        bdxdbDirectMultPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
    }
    return bdxdbDirectMultPtr_()[patchI];
}


const vectorField& objective::boundaryEdgeMultiplier
(
    const label patchI,
    const label edgeI
)
{
    if (bdxdbDirectMultPtr_.empty())
    {
        FatalErrorInFunction
            << "Unallocated boundaryEdgeMultiplier field"
            << exit(FatalError);
    }
    return bEdgeContribution_()[patchI][edgeI];
}


const fvPatchTensorField& objective::boundarydJdStress(const label patchI)
{
    if (bdJdStressPtr_.empty())
    {
        bdJdStressPtr_.reset(createZeroBoundaryPtr<tensor>(mesh_));
    }
    return bdJdStressPtr_()[patchI];
}


const boundaryVectorField& objective::boundarydJdb()
{
    if (bdJdbPtr_.empty())
    {
        bdJdbPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
    }
    return bdJdbPtr_();
}


const boundaryVectorField& objective::dSdbMultiplier()
{
    if (bdSdbMultPtr_.empty())
    {
        bdSdbMultPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
    }
    return bdSdbMultPtr_();
}


const boundaryVectorField& objective::dndbMultiplier()
{
    if (bdndbMultPtr_.empty())
    {
        bdndbMultPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
    }
    return bdndbMultPtr_();
}


const boundaryVectorField& objective::dxdbMultiplier()
{
    if (bdxdbMultPtr_.empty())
    {
        bdxdbMultPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
    }
    return bdxdbMultPtr_();
}


const boundaryVectorField& objective::dxdbDirectMultiplier()
{
    if (bdxdbDirectMultPtr_.empty())
    {
        bdxdbDirectMultPtr_.reset(createZeroBoundaryPtr<vector>(mesh_));
    }
    return bdxdbDirectMultPtr_();
}


const vectorField3& objective::boundaryEdgeMultiplier()
{
    if (bdxdbDirectMultPtr_.empty())
    {
        FatalErrorInFunction
            << "Unallocated boundaryEdgeMultiplier field"
            << endl << endl
            << exit(FatalError);
    }
    return bEdgeContribution_();
}


const boundaryTensorField& objective::boundarydJdStress()
{
    if (bdJdStressPtr_.empty())
    {
        bdJdStressPtr_.reset(createZeroBoundaryPtr<tensor>(mesh_));
    }
    return bdJdStressPtr_();
}


const volScalarField& objective::divDxDbMultiplier()
{
    if (divDxDbMultPtr_.empty())
    {
        // If pointer is not set, set it to a zero field
        divDxDbMultPtr_.reset
        (
            createZeroFieldPtr<scalar>
            (
                mesh_,
                ("divDxDbMult"+objectiveName_),
                // Variable dimensions!!
                // Dummy dimensionless. Only the internalField will be used
                dimless
            )
        );
    }
    return divDxDbMultPtr_();
}


const volTensorField& objective::gradDxDbMultiplier()
{
    if (gradDxDbMultPtr_.empty())
    {
        // If pointer is not set, set it to a zero field
        gradDxDbMultPtr_.reset
        (
            createZeroFieldPtr<tensor>
            (
                mesh_,
                ("gradDxDbMult"+objectiveName_),
                // Variable dimensions!!
                dimensionSet(pow2(dimLength)/pow3(dimTime))
            )
        );
    }
    return gradDxDbMultPtr_();
}


void objective::write() const
{
    if (Pstream::master())
    {
        // File is opened only upon invocation of the write function
        // in order to avoid various instantiations of the same objective
        // opening the same file
        if (objFunctionFilePtr_.empty())
        {
            setObjectiveFilePtr();
        }

        objFunctionFilePtr_() << mesh_.time().value() << tab << J_ << endl;
    }
}


void objective::writeInstantaneousValue() const
{
    if (Pstream::master())
    {
        // File is opened only upon invocation of the write function
        // in order to avoid various instantiations of the same objective
        // opening the same file
        if (instantValueFilePtr_.empty())
        {
            setInstantValueFilePtr();
        }

        instantValueFilePtr_() << mesh_.time().value() << tab << J_ << endl;
    }
}


void objective::writeMeanValue() const
{
    if (Pstream::master())
    {
        // Write mean value if necessary
        if (computeMeanFields_)
        {
            // File is opened only upon invocation of the write function
            // in order to avoid various instantiations of the same objective
            // opening the same file
            if (meanValueFilePtr_.empty())
            {
                setMeanValueFilePtr();
            }

            meanValueFilePtr_()
                << mesh_.time().value() << tab << JMean_ << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
