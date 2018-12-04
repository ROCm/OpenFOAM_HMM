/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2018 OpenCFD Ltd.
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

#include "volRegion.H"
#include "volMesh.H"
#include "cellSet.H"
#include "globalMeshData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(volRegion, 0);
}
}


const Foam::Enum
<
    Foam::functionObjects::volRegion::regionTypes
>
Foam::functionObjects::volRegion::regionTypeNames_
({
    { regionTypes::vrtAll, "all" },
    { regionTypes::vrtCellSet, "cellSet" },
    { regionTypes::vrtCellZone, "cellZone" },
});


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::volRegion::writeFileHeader
(
    const writeFile& wf,
    Ostream& file
) const
{
    wf.writeCommented(file, "Region");
    file<< setw(1) << ':' << setw(1) << ' '
        << regionTypeNames_[regionType_] << " " << regionName_ << endl;
    wf.writeHeaderValue(file, "Cells", nCells_);
    wf.writeHeaderValue(file, "Volume", V_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::volRegion::volRegion
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    volMesh_(mesh),
    regionType_
    (
        regionTypeNames_.lookupOrDefault
        (
            "regionType",
            dict,
            regionTypes::vrtAll
        )
    ),
    regionName_(volMesh_.name()),
    regionID_(-1)
{
    read(dict);

    // Cache integral properties of the region for writeFileHeader
    nCells_ = nCells();
    V_ = V();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::volRegion::read
(
    const dictionary& dict
)
{
    regionID_ = -1;
    cellIds_.clear();

    switch (regionType_)
    {
        case vrtCellSet:
        {
            dict.readEntry("name", regionName_);

            cellIds_ = cellSet(volMesh_, regionName_).sortedToc();

            if (nCells() == 0)
            {
                FatalIOErrorInFunction(dict)
                    << regionTypeNames_[regionType_]
                    << "(" << regionName_ << "):" << nl
                    << "    Region has no cells"
                    << exit(FatalIOError);
            }

            break;
        }

        case vrtCellZone:
        {
            dict.readEntry("name", regionName_);

            regionID_ = volMesh_.cellZones().findZoneID(regionName_);

            if (regionID_ < 0)
            {
                FatalIOErrorInFunction(dict)
                    << "Unknown cell zone name: " << regionName_
                    << ". Valid cell zones    : "
                    << flatOutput(volMesh_.cellZones().names())
                    << exit(FatalIOError);
            }

            if (nCells() == 0)
            {
                FatalIOErrorInFunction(dict)
                    << regionTypeNames_[regionType_]
                    << "(" << regionName_ << "):" << nl
                    << "    Region has no cells"
                    << exit(FatalIOError);
            }

            break;
        }

        case vrtAll:
        {
            regionName_= volMesh_.name();
            break;
        }

        default:
        {
            FatalIOErrorInFunction(dict)
                << "Unknown region type. Valid region types are:"
                << flatOutput(regionTypeNames_.names()) << nl
                << exit(FatalIOError);
        }
    }

    return true;
}


const Foam::labelList& Foam::functionObjects::volRegion::cellIDs() const
{
    switch (regionType_)
    {
        case vrtCellSet:
            return cellIds_;
            break;

        case vrtCellZone:
            return volMesh_.cellZones()[regionID_];
            break;

        default:
            break;
    }

    return labelList::null();
}


Foam::label Foam::functionObjects::volRegion::nCells() const
{
    if (regionType_ == vrtAll)
    {
        return volMesh_.globalData().nTotalCells();
    }
    else
    {
        return returnReduce(cellIDs().size(), sumOp<label>());
    }
}


Foam::scalar Foam::functionObjects::volRegion::V() const
{
    if (regionType_ == vrtAll)
    {
        return gSum(volMesh_.V());
    }
    else
    {
        scalar vol = 0;
        for (const label celli : cellIDs())
        {
            vol += volMesh_.V()[celli];
        }

        return returnReduce(vol, sumOp<scalar>());
    }
}


// ************************************************************************* //
