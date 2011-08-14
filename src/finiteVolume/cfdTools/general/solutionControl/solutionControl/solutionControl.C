/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    const dictionary& solnDict = this->dict();

    // Read solution controls
    nNonOrthCorr_ =
        solnDict.lookupOrDefault<label>("nNonOrthogonalCorrectors", 0);
    momentumPredictor_ = solnDict.lookupOrDefault("momentumPredictor", true);
    transonic_ = solnDict.lookupOrDefault("transonic", false);

    // Read residual information
    const dictionary residualDict(solnDict.subOrEmptyDict("residualControl"));
    DynamicList<fieldData> data(residualDict.toc().size());
    wordHashSet fieldNames;

    forAllConstIter(dictionary, residualDict, iter)
    {
        if (fieldNames.insert(iter().keyword()))
        {
            fieldData fd;
            fd.name = iter().keyword().c_str();

            if (absTolOnly)
            {
                fd.absTol = readScalar(residualDict.lookup(iter().keyword()));
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
                    FatalErrorIn("bool Foam::solutionControl::read()")
                        << "Residual data for " << iter().keyword()
                        << " must be specified as a dictionary";
                }
            }

            data.append(fd);
        }
    }

    residualControl_.transfer(data);
}


Foam::label Foam::solutionControl::applyToField(const word& fieldName) const
{
    forAll(residualControl_, i)
    {
        if (residualControl_[i].name.match(fieldName))
        {
            return i;
        }
    }

    return -1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solutionControl::solutionControl(fvMesh& mesh, const word& algorithmName)
:
    mesh_(mesh),
    residualControl_(),
    algorithmName_(algorithmName),
    nNonOrthCorr_(0),
    momentumPredictor_(true),
    transonic_(false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solutionControl::~solutionControl()
{}


// ************************************************************************* //
