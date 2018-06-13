/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenCFD Ltd.
     \\/     M anipulation  | Copyright (C) 2015 IH-Cantabria
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

#include "waveGenerationModel.H"
#include "mathematicalConstants.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace waveModels
{
    defineTypeNameAndDebug(waveGenerationModel, 0);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::waveModels::waveGenerationModel::readWaveHeight() const
{
    scalar h(readScalar(lookup("waveHeight")));
    if (h < 0)
    {
        FatalIOErrorInFunction(*this)
            << "Wave height must be greater than zero.  Supplied"
            << " value waveHeight = " << h
            << exit(FatalIOError);
    }

    return h;
}


Foam::scalar Foam::waveModels::waveGenerationModel::readWaveAngle() const
{
    scalar angle(readScalar(lookup("waveAngle")));
    return angle* mathematical::pi/180;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveModels::waveGenerationModel::waveGenerationModel
(
    const dictionary& dict,
    const fvMesh& mesh,
    const polyPatch& patch,
    const bool readFields
)
:
    waveModel(dict, mesh, patch, false)
{
    if (readFields)
    {
        readDict(dict);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::waveModels::waveGenerationModel::~waveGenerationModel()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::waveModels::waveGenerationModel::readDict
(
    const dictionary& overrideDict
)
{
    if (waveModel::readDict(overrideDict))
    {
        lookup("activeAbsorption") >> activeAbsorption_;

        return true;
    }

    return false;
}


void Foam::waveModels::waveGenerationModel::info(Ostream& os) const
{
    waveModel::info(os);
}


// ************************************************************************* //
