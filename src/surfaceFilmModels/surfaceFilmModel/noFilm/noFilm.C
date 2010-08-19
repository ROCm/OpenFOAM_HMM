/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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

#include "noFilm.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace surfaceFilmModels
    {
        defineTypeNameAndDebug(noFilm, 0);
        addToRunTimeSelectionTable(surfaceFilmModel, noFilm, mesh);
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::surfaceFilmModels::noFilm::read()
{
    if (surfaceFilmModel::read())
    {
        // no additional info to read
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceFilmModels::noFilm::noFilm
(
    const word& modelType,
    const fvMesh& mesh,
    const dimensionedVector& g
)
:
    surfaceFilmModel(modelType, mesh, g)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceFilmModels::noFilm::~noFilm()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::surfaceFilmModels::noFilm::evolveFilm()
{
    // do nothing
}


const Foam::fvMesh& Foam::surfaceFilmModels::noFilm::film() const
{
    FatalErrorIn("const fvMesh& noFilm::film() const")
        << "Cannot return film for noFilm model" << abort(FatalError);

    return reinterpret_cast<const fvMesh&>(null);
}


const Foam::labelList&
Foam::surfaceFilmModels::noFilm::filmBottomPatchIDs() const
{
    return labelList::null();
}


const Foam::labelList& Foam::surfaceFilmModels::noFilm::filmTopPatchIDs() const
{
    return labelList::null();
}


const Foam::labelList& Foam::surfaceFilmModels::noFilm::primaryPatchIDs() const
{
    return labelList::null();
}


bool Foam::surfaceFilmModels::noFilm::isFilmPatch(const label) const
{
    return false;
}


void Foam::surfaceFilmModels::noFilm::addSources
(
    const label,
    const label,
    const scalar,
    const vector&,
    const scalar,
    const scalar
)
{
    // do nothing
}


const Foam::volVectorField& Foam::surfaceFilmModels::noFilm::U() const
{
    FatalErrorIn
    (
        "const volScalarField& noFilm::U() const"
    )   << "U field not available for " << type() << abort(FatalError);

    return volVectorField::null();
}


const Foam::volScalarField& Foam::surfaceFilmModels::noFilm::rho() const
{
    FatalErrorIn
    (
        "const volScalarField& noFilm::rho() const"
    )   << "rho field not available for " << type() << abort(FatalError);

    return volScalarField::null();
}


const Foam::volScalarField& Foam::surfaceFilmModels::noFilm::T() const
{
    FatalErrorIn
    (
        "const Foam::volScalarField& Foam::noFilm::T() const"
    )   << "T field not available for " << type() << abort(FatalError);

    return volScalarField::null();
}


const Foam::volScalarField& Foam::surfaceFilmModels::noFilm::cp() const
{
    FatalErrorIn
    (
        "const volScalarField& noFilm::cp() const"
    )   << "cp field not available for " << type() << abort(FatalError);

    return volScalarField::null();
}


const Foam::volScalarField&
Foam::surfaceFilmModels::noFilm::massForPrimary() const
{
    FatalErrorIn
    (
        "const volScalarField& noFilm::massForPrimary() const"
    )   << "massForPrimary field not available for " << type()
        << abort(FatalError);

    return volScalarField::null();
}


const Foam::volScalarField&
Foam::surfaceFilmModels::noFilm::diametersForPrimary() const
{
    FatalErrorIn
    (
        "const volScalarField& noFilm::diametersForPrimary() const"
    )   << "diametersForPrimary field not available for " << type()
        << abort(FatalError);

    return volScalarField::null();
}


void Foam::surfaceFilmModels::noFilm::info() const
{
    // do nothing
}


// ************************************************************************* //
