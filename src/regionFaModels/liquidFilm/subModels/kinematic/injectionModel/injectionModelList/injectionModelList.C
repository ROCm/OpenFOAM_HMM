/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "injectionModelList.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace regionModels
{
namespace areaSurfaceFilmModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

injectionModelList::injectionModelList(liquidFilmBase& film)
:
    PtrList<injectionModel>(),
    filmSubModelBase(film)
{}


injectionModelList::injectionModelList
(
    liquidFilmBase& film,
    const dictionary& dict
)
:
    PtrList<injectionModel>(),
    filmSubModelBase
    (
        "injectionModelList",
        film,
        dict,
        "injectionModelList",
        "injectionModelList"
    ),
    massInjected_(Zero)
{
    const wordList activeModels(dict.lookup("injectionModels"));

    wordHashSet models(activeModels);

    Info<< "    Selecting film injection models" << endl;
    if (models.size())
    {
        this->setSize(models.size());

        label i = 0;
        for (const word& model : models)
        {
            set(i, injectionModel::New(film, dict, model));
            i++;
        }
    }
    else
    {
        Info<< "        none" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

injectionModelList::~injectionModelList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void injectionModelList::correct
(
    scalarField& availableMass,
    volScalarField& massToInject,
    volScalarField& diameterToInject
)
{
    const label patchi = film().patchID();

    // Correct models that accumulate mass and diameter transfers
    forAll(*this, i)
    {
        injectionModel& im = operator[](i);
        im.correct
        (
            availableMass,
            massToInject.boundaryFieldRef()[patchi],
            diameterToInject.boundaryFieldRef()[patchi]
        );
    }

    massInjected_ += gSum(massToInject.boundaryField()[patchi]);
}


void injectionModelList::info(Ostream& os)
{
    const polyBoundaryMesh& pbm = film().primaryMesh().boundaryMesh();

    scalar injectedMass = 0;
    scalar patchInjectedMasses = 0;

    forAll(*this, i)
    {
        const injectionModel& im = operator[](i);
        injectedMass += im.injectedMassTotal();
        im.patchInjectedMassTotals(patchInjectedMasses);
    }

    os  << indent << "injected mass      = " << injectedMass << nl;

    const label patchi = film().patchID();

    if (mag(patchInjectedMasses) > VSMALL)
    {
        os  << indent << indent << "from patch " << pbm[patchi].name()
            << " = " << patchInjectedMasses << nl;
    }

    scalar mass0(Zero);
    this->getBaseProperty("massInjected", mass0);

    scalar mass(massInjected_);
    mass += mass0;

    Info<< indent << "  - patch: " << pbm[patchi].name() << "  "
        << mass << endl;

    if (film().primaryMesh().time().writeTime())
    {
        setBaseProperty("massInjected", mass);
        massInjected_ = 0.0;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace areaSurfaceFilmModels
} // End namespace regionModels
} // End namespace Foam

// ************************************************************************* //
