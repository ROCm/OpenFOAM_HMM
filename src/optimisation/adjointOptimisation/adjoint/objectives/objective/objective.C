/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
    Copyright (C) 2019 OpenCFD Ltd.
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
    nullified_(false),

    J_(Zero),
    JMean_(Zero),
    weight_(Zero),

    integrationStartTimePtr_(nullptr),
    integrationEndTimePtr_(nullptr),

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
    // Read integration start and end times, if present.
    // For unsteady runs only
    if (dict.found("integrationStartTime"))
    {
        integrationStartTimePtr_.reset
        (
            new scalar(dict.get<scalar>("integrationStartTime"))
        );
    }
    if (dict.found("integrationEndTime"))
    {
        integrationEndTimePtr_.reset
        (
            new scalar(dict.get<scalar>("integrationEndTime"))
        );
    }

    // Read JMean from dictionary, if present
    IOobject headObjectiveIODict
    (
        "objectiveDict" + objectiveName_,
        mesh_.time().timeName(),
        "uniform",
        mesh_,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    );
    if (headObjectiveIODict.typeHeaderOk<IOdictionary>(false))
    {
        JMean_ = IOdictionary(headObjectiveIODict).get<scalar>("JMean");
    }
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
        FatalIOErrorInLookup
        (
            dict,
            "objective",
            objectiveType,
            *objectiveConstructorTablePtr_
        ) << exit(FatalIOError);
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


scalar objective::JCycle() const
{
    scalar J(J_);
    if
    (
        computeMeanFields_
     || (hasIntegrationStartTime() && hasIntegrationEndTime())
    )
    {
        J = JMean_;
    }
    return J;
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


void objective::accumulateJMean()
{
    if (hasIntegrationStartTime() && hasIntegrationEndTime())
    {
        const scalar time = mesh_.time().value();
        if (isWithinIntegrationTime())
        {
            const scalar dt = mesh_.time().deltaT().value();
            const scalar elapsedTime = time - integrationStartTimePtr_();
            const scalar denom = elapsedTime + dt;
            JMean_ = (JMean_*elapsedTime + J_*dt)/denom;
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unallocated integration start or end time"
            << exit(FatalError);
    }
}


scalar objective::weight() const
{
    return weight_;
}


bool objective::isWithinIntegrationTime() const
{
    if (hasIntegrationStartTime() && hasIntegrationEndTime())
    {
        const scalar time = mesh_.time().value();
        return
            (
                time >= integrationStartTimePtr_()
             && time <= integrationEndTimePtr_()
            );
    }
    else
    {
        FatalErrorInFunction
            << "Unallocated integration start or end time"
            << exit(FatalError);
    }
    return false;
}


void objective::incrementIntegrationTimes(const scalar timeSpan)
{
    if (hasIntegrationStartTime() && hasIntegrationEndTime())
    {
        integrationStartTimePtr_() += timeSpan;
        integrationEndTimePtr_() += timeSpan;
    }
    else
    {
        FatalErrorInFunction
            << "Unallocated integration start or end time"
            << exit(FatalError);
    }
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


void objective::nullify()
{
    if (!nullified_)
    {
        if (hasdJdb())
        {
            dJdbPtr_() == dimensionedScalar(dJdbPtr_().dimensions(), Zero);
        }
        if (hasBoundarydJdb())
        {
            bdJdbPtr_() == vector::zero;
        }
        if (hasdSdbMult())
        {
            bdSdbMultPtr_() == vector::zero;
        }
        if (hasdndbMult())
        {
            bdndbMultPtr_() == vector::zero;
        }
        if (hasdxdbMult())
        {
            bdxdbMultPtr_() == vector::zero;
        }
        if (hasdxdbDirectMult())
        {
            bdxdbDirectMultPtr_() == vector::zero;
        }
        if (hasBoundaryEdgeContribution())
        {
            for (Field<vectorField>& field : bEdgeContribution_())
            {
                field = vector::zero;
            }
        }
        if (hasBoundarydJdStress())
        {
            bdJdStressPtr_() == tensor::zero;
        }
        if (hasDivDxDbMult())
        {
            divDxDbMultPtr_() ==
                dimensionedScalar(divDxDbMultPtr_().dimensions(), Zero);
        }
        if (hasGradDxDbMult())
        {
            gradDxDbMultPtr_() ==
                dimensionedTensor(gradDxDbMultPtr_().dimensions(), Zero);
        }

        nullified_ = true;
    }
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
        // Covers both steady and unsteady runs
        if
        (
            computeMeanFields_
         || (hasIntegrationStartTime() && hasIntegrationEndTime())
        )
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
    // Write mean value under time/uniform, to allow for lineSearch to work
    // appropriately in continuation runs, when field averaging is used
    IOdictionary objectiveDict
    (
        IOobject
        (
            "objectiveDict" + objectiveName_,
            mesh_.time().timeName(),
            "uniform",
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        )
    );
    objectiveDict.add<scalar>("JMean", JMean_);
    objectiveDict.regIOobject::write();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
