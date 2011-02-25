/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2011 OpenCFD Ltd.
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

#include "cloudSolution.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cloudSolution::cloudSolution
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    dict_(dict),
    active_(dict.lookup("active")),
    transient_(false),
    calcFrequency_(1),
    maxCo_(0.3),
    iter_(1),
    deltaT_(0.0),
    coupled_(false),
    cellValueSourceCorrection_(false),
    maxTrackTime_(0.0),
    resetSourcesOnStartup_(true)
{
    if (active_)
    {
        read();
    }
}


Foam::cloudSolution::cloudSolution
(
    const cloudSolution& cs
)
:
    mesh_(cs.mesh_),
    dict_(cs.dict_),
    active_(cs.active_),
    transient_(cs.transient_),
    calcFrequency_(cs.calcFrequency_),
    maxCo_(cs.maxCo_),
    iter_(cs.iter_),
    deltaT_(cs.deltaT_),
    coupled_(cs.coupled_),
    cellValueSourceCorrection_(cs.cellValueSourceCorrection_),
    maxTrackTime_(cs.maxTrackTime_),
    resetSourcesOnStartup_(cs.resetSourcesOnStartup_)
{}


Foam::cloudSolution::cloudSolution
(
    const fvMesh& mesh
)
:
    mesh_(mesh),
    dict_(dictionary::null),
    active_(false),
    transient_(false),
    calcFrequency_(0),
    maxCo_(GREAT),
    iter_(0),
    deltaT_(0.0),
    coupled_(false),
    cellValueSourceCorrection_(false),
    maxTrackTime_(0.0),
    resetSourcesOnStartup_(false)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cloudSolution::~cloudSolution()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cloudSolution::read()
{
    dict_.lookup("transient") >> transient_;
    dict_.lookup("coupled") >> coupled_;
    dict_.lookup("cellValueSourceCorrection") >> cellValueSourceCorrection_;

    if (steadyState())
    {
        dict_.lookup("calcFrequency") >> calcFrequency_;
        dict_.lookup("maxCo") >> maxCo_;
        dict_.lookup("maxTrackTime") >> maxTrackTime_;
        dict_.subDict("sourceTerms").lookup("resetOnStartup")
            >> resetSourcesOnStartup_;
    }
}


Foam::scalar Foam::cloudSolution::relaxCoeff(const word& fieldName) const
{
    return readScalar(sourceTermDict().subDict(fieldName).lookup("alpha"));
}


bool Foam::cloudSolution::semiImplicit(const word& fieldName) const
{
    return readBool(sourceTermDict().subDict(fieldName).lookup("semiImplicit"));
}


bool Foam::cloudSolution::solveThisStep() const
{
    return
        active_
     && (
            mesh_.time().outputTime()
         || (mesh_.time().timeIndex() % calcFrequency_ == 0)
        );
}


bool Foam::cloudSolution::canEvolve()
{
    // Set the calculation time step
    if (transient_)
    {
        deltaT_ = mesh_.time().deltaTValue();
    }
    else
    {
        deltaT_ = maxTrackTime_;
    }

    return solveThisStep();
}


bool Foam::cloudSolution::output() const
{
    return active_ && mesh_.time().outputTime();
}


// ************************************************************************* //
