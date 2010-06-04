/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
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

#include "residualControl.H"
#include "fvMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::residualControl, 0);


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::residualControl::checkCriteria(const bool verbose) const
{
    bool achieved = true;
    const fvMesh& mesh = static_cast<const fvMesh&>(obr_);
    const dictionary& solverDict = mesh.solverPerformanceDict();

    forAllConstIter(dictionary, solverDict, iter)
    {
        const word& variableName = iter().keyword();
        scalar maxResidual;

        if (maxResiduals_.readIfPresent(variableName, maxResidual))
        {
            // use the residual from the first solution
            const scalar eqnResidual =
                List<lduMatrix::solverPerformance>
                (
                    iter().stream()
                ).first().initialResidual();

            achieved = achieved && (eqnResidual < maxResidual);

            if (verbose)
            {
                Info<< "    " << variableName
                    << ": requested max residual = " << maxResidual
                    << ", eqn residual = " << eqnResidual << nl;
            }
        }
    }

    return achieved;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::residualControl::residualControl
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool
)
:
    name_(name),
    obr_(obr),
    active_(true),
    maxResiduals_(),
    criteriaSatisfied_(false)
{
    // Only active if a fvMesh is available
    if (isA<fvMesh>(obr_))
    {
        read(dict);
    }
    else
    {
        active_ = false;
        WarningIn
        (
            "residualControl::residualControl"
            "("
                "const word&, "
                "const objectRegistry&, "
                "const dictionary&, "
                "const bool"
            ")"
        )   << "No fvMesh available, deactivating."
            << nl << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::residualControl::~residualControl()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::residualControl::read(const dictionary& dict)
{
    if (active_)
    {
        maxResiduals_ = dict.subDict("maxResiduals");
    }
}


void Foam::residualControl::execute()
{
    if (active_)
    {
        criteriaSatisfied_ = checkCriteria(false);

        if (criteriaSatisfied_)
        {
            Info<< "Convergence criteria satisfied - finalising run" << nl
                << endl;

            checkCriteria(true);

            Info<< endl;

            const fvMesh& mesh = static_cast<const fvMesh&>(obr_);
            Time& time = const_cast<Time&>(mesh.time());
            time.writeAndEnd();
        }
    }
}


void Foam::residualControl::end()
{
    if (active_)
    {
        if (criteriaSatisfied_)
        {
            Info<< "Residual control criteria satisfied" << nl;
        }
        else
        {
            Info<< "Residual control criteria not satisfied" << nl;
        }
    }
}


void Foam::residualControl::write()
{
    // do nothing
}


// ************************************************************************* //
