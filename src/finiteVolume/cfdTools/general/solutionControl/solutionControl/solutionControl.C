/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "solutionControl.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(solutionControl, 0);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::solutionControl::read(const bool absTolOnly)
{
    const dictionary& solutionDict = this->dict();

    // Read solution controls
    nNonOrthCorr_ =
        solutionDict.lookupOrDefault<label>("nNonOrthogonalCorrectors", 0);
    momentumPredictor_ =
        solutionDict.lookupOrDefault("momentumPredictor", true);
    transonic_ = solutionDict.lookupOrDefault("transonic", false);
    consistent_ = solutionDict.lookupOrDefault("consistent", false);
    frozenFlow_ = solutionDict.lookupOrDefault("frozenFlow", false);

    // Read residual information
    const dictionary residualDict
    (
        solutionDict.subOrEmptyDict("residualControl")
    );

    DynamicList<fieldData> data(residualControl_);

    forAllConstIter(dictionary, residualDict, iter)
    {
        const word& fName = iter().keyword();
        const label fieldi = applyToField(fName, false);
        if (fieldi == -1)
        {
            fieldData fd;
            fd.name = fName.c_str();

            if (absTolOnly)
            {
                fd.absTol = readScalar(residualDict.lookup(fName));
                fd.relTol = -1;
                fd.initialResidual = -1;
            }
            else
            {
                if (iter().isDict())
                {
                    const dictionary& fieldDict(iter().dict());
                    fd.absTol = readScalar(fieldDict.lookup("tolerance"));
                    fd.relTol = readScalar(fieldDict.lookup("relTol"));
                    fd.initialResidual = 0.0;
                }
                else
                {
                    FatalErrorInFunction
                        << "Residual data for " << iter().keyword()
                        << " must be specified as a dictionary"
                        << exit(FatalError);
                }
            }

            data.append(fd);
        }
        else
        {
            fieldData& fd = data[fieldi];
            if (absTolOnly)
            {
                fd.absTol = readScalar(residualDict.lookup(fName));
            }
            else
            {
                if (iter().isDict())
                {
                    const dictionary& fieldDict(iter().dict());
                    fd.absTol = readScalar(fieldDict.lookup("tolerance"));
                    fd.relTol = readScalar(fieldDict.lookup("relTol"));
                }
                else
                {
                    FatalErrorInFunction
                        << "Residual data for " << iter().keyword()
                        << " must be specified as a dictionary"
                        << exit(FatalError);
                }
            }
        }
    }

    residualControl_.transfer(data);

    if (debug)
    {
        forAll(residualControl_, i)
        {
            const fieldData& fd = residualControl_[i];
            Info<< "residualControl[" << i << "]:" << nl
                << "    name     : " << fd.name << nl
                << "    absTol   : " << fd.absTol << nl
                << "    relTol   : " << fd.relTol << nl
                << "    iniResid : " << fd.initialResidual << endl;
        }
    }
}


void Foam::solutionControl::read()
{
    read(false);
}


Foam::label Foam::solutionControl::applyToField
(
    const word& fieldName,
    const bool useRegEx
) const
{
    forAll(residualControl_, i)
    {
        if (residualControl_[i].name.match(fieldName, !useRegEx))
        {
            return i;
        }
    }

    return -1;
}


void Foam::solutionControl::storePrevIterFields() const
{
//    storePrevIter<label>();
    storePrevIter<scalar>();
    storePrevIter<vector>();
    storePrevIter<sphericalTensor>();
    storePrevIter<symmTensor>();
    storePrevIter<tensor>();
}


Foam::Pair<Foam::scalar> Foam::solutionControl::maxResidual
(
    const entry& solverPerfDictEntry
) const
{
    return maxResidual(mesh_, solverPerfDictEntry);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Type>
bool Foam::solutionControl::maxTypeResidual
(
    const fvMesh& fvmesh,
    const entry& solverPerfDictEntry,
    Pair<scalar>& residuals
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    const word& fieldName = solverPerfDictEntry.keyword();

    if (fvmesh.foundObject<fieldType>(fieldName))
    {
        const List<SolverPerformance<Type>> sp(solverPerfDictEntry.stream());

        residuals.first() = cmptMax(sp.first().initialResidual());
        residuals.last()  = cmptMax(sp.last().initialResidual());

        return true;
    }

    return false;
}


Foam::Pair<Foam::scalar> Foam::solutionControl::maxResidual
(
    const fvMesh& fvmesh,
    const entry& solverPerfDictEntry
)
{
    Pair<scalar> residuals(0,0);

    // Check with builtin short-circuit
    const bool ok =
    (
        maxTypeResidual<scalar>(fvmesh, solverPerfDictEntry, residuals)
     || maxTypeResidual<vector>(fvmesh, solverPerfDictEntry, residuals)
     || maxTypeResidual<sphericalTensor>(fvmesh, solverPerfDictEntry, residuals)
     || maxTypeResidual<symmTensor>(fvmesh, solverPerfDictEntry, residuals)
     || maxTypeResidual<tensor>(fvmesh, solverPerfDictEntry, residuals)
    );

    if (!ok && solutionControl::debug)
    {
        Info<<"no residual for " << solverPerfDictEntry.keyword()
            << " on mesh " << fvmesh.name() << nl;
    }

    return residuals;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solutionControl::solutionControl(fvMesh& mesh, const word& algorithmName)
:
    IOobject
    (
        "solutionControl",
        mesh.time().timeName(),
        mesh
    ),
    mesh_(mesh),
    residualControl_(),
    algorithmName_(algorithmName),
    nNonOrthCorr_(0),
    momentumPredictor_(true),
    transonic_(false),
    consistent_(false),
    frozenFlow_(false),
    corr_(0),
    corrNonOrtho_(0)
{}


// ************************************************************************* //
