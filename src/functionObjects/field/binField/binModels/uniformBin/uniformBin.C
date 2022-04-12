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

#include "uniformBin.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace binModels
{
    defineTypeNameAndDebug(uniformBin, 0);
    addToRunTimeSelectionTable(binModel, uniformBin, dictionary);
}
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::binModels::uniformBin::initialise()
{
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    // Determine extents of patches in a given coordinate system
    vector geomMin(GREAT, GREAT, GREAT);
    vector geomMax(-GREAT, -GREAT, -GREAT);

    for (const label patchi : patchSet_)
    {
        const polyPatch& pp = pbm[patchi];
        const vectorField ppcs(coordSysPtr_->localPosition(pp.faceCentres()));

        for (direction i = 0; i < vector::nComponents; ++i)
        {
            geomMin[i] = min(min(ppcs.component(i)), geomMin[i]);
            geomMax[i] = max(max(ppcs.component(i)), geomMax[i]);
        }
    }

    for (const label zonei : cellZoneIDs_)
    {
        const cellZone& cZone = mesh_.cellZones()[zonei];
        const vectorField d
        (
            coordSysPtr_->localPosition(vectorField(mesh_.C(), cZone))
        );

        for (direction i = 0; i < vector::nComponents; ++i)
        {
            geomMin[i] = min(min(d.component(i)), geomMin[i]);
            geomMax[i] = max(max(d.component(i)), geomMax[i]);
        }
    }

    reduce(geomMin, minOp<vector>());
    reduce(geomMax, maxOp<vector>());

    for (direction i = 0; i < vector::nComponents; ++i)
    {
        // Slightly boost max so that region of interest is fully within bounds
        geomMax[i] = 1.0001*(geomMax[i] - geomMin[i]) + geomMin[i];

        // Use geometry limits if not specified by the user
        if (binMinMax_[i][0] == GREAT) binMinMax_[i][0] = geomMin[i];
        if (binMinMax_[i][1] == GREAT) binMinMax_[i][1] = geomMax[i];

        if (binMinMax_[i][0] > binMinMax_[i][1])
        {
            FatalErrorInFunction
                << "Max bounds must be greater than min bounds" << nl
                << "    direction   = " << i << nl
                << "    min         = " << binMinMax_[i][0] << nl
                << "    max         = " << binMinMax_[i][1] << nl
                << exit(FatalError);
        }

        //- Compute bin widths in binning directions
        binW_[i] = (binMinMax_[i][1] - binMinMax_[i][0])/scalar(nBins_[i]);

        if (binW_[i] <= 0)
        {
            FatalErrorInFunction
                << "Bin widths must be greater than zero" << nl
                << "    direction = " << i << nl
                << "    min bound = " << binMinMax_[i][0] << nl
                << "    max bound = " << binMinMax_[i][1] << nl
                << "    bin width = " << binW_[i]
                << exit(FatalError);
        }
    }

    setBinsAddressing();
}


Foam::labelList Foam::binModels::uniformBin::binAddr(const vectorField& d) const
{
    labelList binIndices(d.size(), -1);

    forAll(d, i)
    {
        // Avoid elements outside of the bin
        bool faceInside = true;
        for (direction j = 0; j < vector::nComponents; ++j)
        {
            if (d[i][j] < binMinMax_[j][0] || d[i][j] > binMinMax_[j][1])
            {
                faceInside = false;
                break;
            }
        }

        if (faceInside)
        {
            // Find the bin division corresponding to the element
            Vector<label> n(Zero);
            for (direction j = 0; j < vector::nComponents; ++j)
            {
                n[j] = floor((d[i][j] - binMinMax_[j][0])/binW_[j]);
                n[j] = min(max(n[j], 0), nBins_[j] - 1);
            }

            // Order: (e1, e2, e3), the first varies the fastest
            binIndices[i] = n[0] + nBins_[0]*n[1] + nBins_[0]*nBins_[1]*n[2];
        }
        else
        {
            binIndices[i] = -1;
        }
    }

    return binIndices;
}


void Foam::binModels::uniformBin::setBinsAddressing()
{
    faceToBin_.setSize(mesh_.nBoundaryFaces());
    faceToBin_ = -1;

    forAllIters(patchSet_, iter)
    {
        const polyPatch& pp = mesh_.boundaryMesh()[iter()];
        const label i0 = pp.start() - mesh_.nInternalFaces();

        SubList<label>(faceToBin_, pp.size(), i0) =
            binAddr(coordSysPtr_->localPosition(pp.faceCentres()));
    }

    cellToBin_.setSize(mesh_.nCells());
    cellToBin_ = -1;

    for (const label zonei : cellZoneIDs_)
    {
        const cellZone& cZone = mesh_.cellZones()[zonei];
        labelList bins
        (
            binAddr(coordSysPtr_->localPosition(vectorField(mesh_.C(), cZone)))
        );

        forAll(cZone, i)
        {
            const label celli = cZone[i];
            cellToBin_[celli] = bins[i];
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::binModels::uniformBin::uniformBin
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& outputPrefix
)
:
    binModel(dict, mesh, outputPrefix),
    nBins_(Zero),
    binW_(Zero),
    binMinMax_
    (
        vector2D(GREAT, GREAT),
        vector2D(GREAT, GREAT),
        vector2D(GREAT, GREAT)
    )
{
    read(dict);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::binModels::uniformBin::read(const dictionary& dict)
{
    if (!binModel::read(dict))
    {
        return false;
    }

    Info<< "    Activating a set of uniform bins" << endl;

    const dictionary& binDict = dict.subDict("binData");

    nBins_ = binDict.get<Vector<label>>("nBin");

    for (const label n : nBins_)
    {
        nBin_ *= n;
    }

    if (nBin_ <= 0)
    {
        FatalIOErrorInFunction(binDict)
            << "Number of bins must be greater than zero" << nl
            << "    e1 bins = " << nBins_[0] << nl
            << "    e2 bins = " << nBins_[1] << nl
            << "    e3 bins = " << nBins_[2]
            << exit(FatalIOError);
    }

    Info<< "    - Employing:" << nl
        << "        " << nBins_[0] << " e1 bins," << nl
        << "        " << nBins_[1] << " e2 bins," << nl
        << "        " << nBins_[2] << " e3 bins"
        << endl;

    cumulative_ = binDict.getOrDefault<bool>("cumulative", false);
    Info<< "    - cumulative    : " << cumulative_ << endl;
    Info<< "    - decomposePatchValues    : " << decomposePatchValues_ << endl;

    if (binDict.found("minMax"))
    {
        const dictionary& minMaxDict = binDict.subDict("minMax");

        for (direction i = 0; i < vector::nComponents; ++i)
        {
            const word ei("e" + Foam::name(i));

            if (minMaxDict.readIfPresent(ei, binMinMax_[i]))
            {
                Info<< "    - " << ei << " min        : "
                    << binMinMax_[i][0] << nl
                    << "    - " << ei << " max        : "
                    << binMinMax_[i][1] << endl;
            }
        }
    }
    Info<< endl;

    initialise();

    return true;
}


void Foam::binModels::uniformBin::apply()
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
                << endl;
        }
    }

    writtenHeader_ = true;
}


void Foam::binModels::uniformBin::updateMesh(const mapPolyMesh& mpm)
{}


void Foam::binModels::uniformBin::movePoints(const polyMesh& mesh)
{
    setBinsAddressing();
}


// ************************************************************************* //
