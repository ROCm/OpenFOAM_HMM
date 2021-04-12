/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

#include "DMD.H"
#include "DMDModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(DMD, 0);
    addToRunTimeSelectionTable(functionObject, DMD, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::DMD::snapshot()
{
    bool processed = false;
    processed = processed || getSnapshot<scalar>();
    processed = processed || getSnapshot<vector>();
    processed = processed || getSnapshot<sphericalTensor>();
    processed = processed || getSnapshot<symmTensor>();
    processed = processed || getSnapshot<tensor>();

    if (!processed)
    {
        FatalErrorInFunction
            << "    functionObjects::" << type() << " " << name() << ":"
            << " cannot find required input field during snapshot loading: "
            << fieldName_ << nl
            << "    Do you execute required functionObjects"
            << " before executing DMD, e.g. mapFields?"
            << exit(FatalError);
    }
}


void Foam::functionObjects::DMD::initialise()
{
    const label nComps = DMDModelPtr_->nComponents(fieldName_);

    if (patch_.empty())
    {
        nSnap_ = nComps*mesh_.nCells();
    }
    else
    {
        const label patchi = mesh_.boundaryMesh().findPatchID(patch_);

        if (patchi < 0)
        {
            FatalErrorInFunction
                << "Cannot find patch " << patch_
                << exit(FatalError);
        }

        nSnap_ = nComps*(mesh_.C().boundaryField()[patchi]).size();
    }

    const label nSnapTotal = returnReduce(nSnap_, sumOp<label>());

    if (nSnapTotal <= 0)
    {
        FatalErrorInFunction
            << "  # Zero-size input field = " << fieldName_ << " #"
            << exit(FatalError);
    }

    if (nSnap_ > 0)
    {
        z_ = RMatrix(2*nSnap_, 1, Zero);
    }
    else
    {
        z_ = RMatrix(1, 1, Zero);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::DMD::DMD
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    DMDModelPtr_(DMDModel::New(mesh_, name, dict)),
    z_(),
    fieldName_(dict.get<word>("field")),
    patch_(dict.getOrDefault<word>("patch", word::null)),
    nSnap_(0),
    step_(0)
{
    // Check if computation time-step is varying
    if (runTime.isAdjustTimeStep())
    {
        WarningInFunction
            << "  # DMD: Available only for fixed time-step computations. #"
            << endl;
    }

    // Check if mesh topology is changing
    if (mesh_.topoChanging())
    {
        FatalErrorInFunction
            << "  # DMD: Available only for non-changing mesh topology. #"
            << exit(FatalError);
    }

    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::DMD::read(const dictionary& dict)
{
    Info<< type() << " " << name() << ":" << endl;

    if (fvMeshFunctionObject::read(dict) && DMDModelPtr_->read(dict))
    {
        return true;
    }

    return false;
}


bool Foam::functionObjects::DMD::execute()
{
    Log << type() << " " << name() << " execute:" << endl;

    snapshot();

    if (step_ == 1)
    {
        DMDModelPtr_->initialise(z_);
    }

    if (step_ > 1)
    {
        DMDModelPtr_->update(z_);
    }

    ++step_;

    Log << tab << "Execution index = " << step_ << " for field: " << fieldName_
        << endl;

    return true;
}


bool Foam::functionObjects::DMD::write()
{
    if (postProcess)
    {
        return true;
    }

    return end();
}


bool Foam::functionObjects::DMD::end()
{
    if (step_ == 0)
    {
        // Avoid double execution of write() when writeControl=onEnd
        return true;
    }

    Log << type() << " " << name() << " write:" << endl;

    if (step_ < 2)
    {
        WarningInFunction
            << "  # DMD needs at least three snapshots to produce output #"
            << nl
            << "  # Only " << step_ + 1 << " snapshots are available #"
            << nl
            << "  # Skipping DMD output calculation and write #"
            << endl;

        return false;
    }

    z_.clear();

    DMDModelPtr_->fit();

    mesh_.time().printExecutionTime(Info);

    // Restart the incremental orthonormal basis update
    step_ = 0;

    return true;
}


// ************************************************************************* //
