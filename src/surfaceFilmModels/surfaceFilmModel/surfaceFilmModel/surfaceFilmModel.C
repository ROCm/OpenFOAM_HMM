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

#include "surfaceFilmModel.H"
#include "fvc.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace surfaceFilmModels
    {
        defineTypeNameAndDebug(surfaceFilmModel, 0);
        defineRunTimeSelectionTable(surfaceFilmModel, mesh);
    }
}


template<>
const char*
Foam::NamedEnum<Foam::surfaceFilmModels::surfaceFilmModel::thermoModelType, 2>::
names[] =
{
    "constant",
    "singleComponent"
};


const
Foam::NamedEnum<Foam::surfaceFilmModels::surfaceFilmModel::thermoModelType, 2>
    Foam::surfaceFilmModels::surfaceFilmModel::thermoModelTypeNames_;


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::surfaceFilmModels::surfaceFilmModel::read()
{
    if (regIOobject::read())
    {
        active_.readIfPresent("active", *this);
        if (const dictionary* dictPtr = subDictPtr(type() + "Coeffs"))
        {
            coeffs_ <<= *dictPtr;
        }

        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceFilmModels::surfaceFilmModel::surfaceFilmModel
(
    const fvMesh& mesh,
    const dimensionedVector& g
)
:
    IOdictionary
    (
        IOobject
        (
            "surfaceFilmProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    time_(mesh.time()),
    active_(false),
    g_(g),
    filmRegionName_("none"),
    coeffs_(dictionary::null),
    thermoModel_(tmConstant)
{}


Foam::surfaceFilmModels::surfaceFilmModel::surfaceFilmModel
(
    const word& type,
    const fvMesh& mesh,
    const dimensionedVector& g
)
:
    IOdictionary
    (
        IOobject
        (
            "surfaceFilmProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    mesh_(mesh),
    time_(mesh.time()),

    active_(lookup("active")),
    g_(g),
    filmRegionName_(lookup("filmRegionName")),
    coeffs_(subDict(type + "Coeffs")),
    thermoModel_(thermoModelTypeNames_.read(coeffs_.lookup("thermoModel")))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceFilmModels::surfaceFilmModel::~surfaceFilmModel()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

const Foam::fvMesh& Foam::surfaceFilmModels::surfaceFilmModel::mesh() const
{
    return mesh_;
}


const Foam::Time& Foam::surfaceFilmModels::surfaceFilmModel::time() const
{
    return time_;
}


void Foam::surfaceFilmModels::surfaceFilmModel::evolve()
{
    if (active_)
    {
        if (mesh_.changing())
        {
            FatalErrorIn("surfaceFilmModel::evolveFilm()")
                << "Currently not possible to apply surface film model to "
                << "moving mesh cases" << nl << abort(FatalError);
        }

        Info<< "\nEvolving surface film for region " << filmRegionName_
            << endl;

        // Update any input information
        read();

        // Pre-evolve
        preEvolveFilm();

        // Increment the film equations up to the new time level
        evolveFilm();

        // Provide some feedback
        Info<< incrIndent;
        info();
        Info<< endl << decrIndent;
    }
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >
Foam::surfaceFilmModels::surfaceFilmModel::Srho() const
{
    notImplemented
    (
        "Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> > "
        "Foam::surfaceFilmModels::surfaceFilmModel::Srho() const"
    )
    return tmp<DimensionedField<scalar, volMesh> >(NULL);
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >
Foam::surfaceFilmModels::surfaceFilmModel::Srho(const label) const
{
    notImplemented
    (
        "Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> > "
        "Foam::surfaceFilmModels::surfaceFilmModel::Srho(const label) const"
    )
    return tmp<DimensionedField<scalar, volMesh> >(NULL);
}


Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> >
Foam::surfaceFilmModels::surfaceFilmModel::Sh() const
{
    notImplemented
    (
        "Foam::tmp<Foam::DimensionedField<Foam::scalar, Foam::volMesh> > "
        "Foam::surfaceFilmModels::surfaceFilmModel::Sh() const"
    )
    return tmp<DimensionedField<scalar, volMesh> >(NULL);
}


// ************************************************************************* //
