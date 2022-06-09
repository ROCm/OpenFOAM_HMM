/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2022 OpenCFD Ltd.
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

#include "singleDirectionUniformBin.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace binModels
{
    defineTypeNameAndDebug(singleDirectionUniformBin, 0);
    addToRunTimeSelectionTable(binModel, singleDirectionUniformBin, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::binModels::singleDirectionUniformBin::initialise()
{
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    // Determine extents of patches in a given direction
    scalar geomMin = GREAT;
    scalar geomMax = -GREAT;
    for (const label patchi : patchSet_)
    {
        const polyPatch& pp = pbm[patchi];
        const scalarField d(pp.faceCentres() & binDir_);
        geomMin = min(min(d), geomMin);
        geomMax = max(max(d), geomMax);
    }

    for (const label zonei : cellZoneIDs_)
    {
        const cellZone& cZone = mesh_.cellZones()[zonei];
        const vectorField cz(mesh_.C(), cZone);
        const scalarField d(cz & binDir_);

        geomMin = min(min(d), geomMin);
        geomMax = max(max(d), geomMax);
    }

    reduce(geomMin, minOp<scalar>());
    reduce(geomMax, maxOp<scalar>());

    // Slightly boost max so that region of interest is fully within bounds
    geomMax = 1.0001*(geomMax - geomMin) + geomMin;

    // Use geometry limits if not specified by the user
    if (binMin_ == GREAT) binMin_ = geomMin;
    if (binMax_ == GREAT) binMax_ = geomMax;

    binDx_ = (binMax_ - binMin_)/scalar(nBin_);

    if (binDx_ <= 0)
    {
        FatalErrorInFunction
            << "Max bound must be greater than min bound" << nl
            << "    d           = " << binDx_ << nl
            << "    min         = " << binMin_ << nl
            << "    max         = " << binMax_ << nl
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::binModels::singleDirectionUniformBin::singleDirectionUniformBin
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& outputPrefix
)
:
    binModel(dict, mesh, outputPrefix),
    binDx_(0),
    binMin_(GREAT),
    binMax_(GREAT),
    binDir_(Zero)
{
    read(dict);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::binModels::singleDirectionUniformBin::read(const dictionary& dict)
{
    if (!binModel::read(dict))
    {
        return false;
    }

    Info<< "    Activating a set of single-direction bins" << endl;

    const dictionary& binDict = dict.subDict("binData");

    nBin_ = binDict.getCheck<label>("nBin", labelMinMax::ge(1));

    Info<< "    Employing " << nBin_ << " bins" << endl;
    if (binDict.readIfPresent("min", binMin_))
    {
        Info<< "    - min        : " << binMin_ << endl;
    }
    if (binDict.readIfPresent("max", binMax_))
    {
        Info<< "    - max        : " << binMax_ << endl;
    }

    cumulative_ = binDict.getOrDefault<bool>("cumulative", false);
    Info<< "    - cumulative    : " << cumulative_ << endl;
    Info<< "    - decomposePatchValues    : " << decomposePatchValues_ << endl;

    binDir_ = binDict.get<vector>("direction");
    binDir_.normalise();

    if (mag(binDir_) == 0)
    {
        FatalIOErrorInFunction(dict)
            << "Input direction should not be zero valued" << nl
            << "    direction = " << binDir_ << nl
            << exit(FatalIOError);
    }

    Info<< "    - direction     : " << binDir_ << nl << endl;

    initialise();

    return true;
}


void Foam::binModels::singleDirectionUniformBin::apply()
{
    forAll(fieldNames_, i)
    {
        const bool ok =
            processField<scalar>(i)
         || processField<vector>(i)
         || processField<sphericalTensor>(i)
         || processField<symmTensor>(i)
         || processField<tensor>(i);

        if (!ok)
        {
            WarningInFunction
                << "Unable to find field " << fieldNames_[i]
                << ". Avaliable objects are "
                << mesh_.objectRegistry::sortedToc()
                << endl;
        }
    }

    writtenHeader_ = true;
}


// ************************************************************************* //
