/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

#include "cellSetOption.H"
#include "cellSet.H"
#include "cellBitSet.H"
#include "volFields.H"
#include "cellCellStencilObject.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace fv
    {
        defineTypeNameAndDebug(cellSetOption, 0);
    }
}


const Foam::Enum
<
    Foam::fv::cellSetOption::selectionModeType
>
Foam::fv::cellSetOption::selectionModeTypeNames_
({
    { selectionModeType::smAll, "all" },
    { selectionModeType::smGeometric, "geometric" },
    { selectionModeType::smPoints, "points" },
    { selectionModeType::smCellSet, "cellSet" },
    { selectionModeType::smCellZone, "cellZone" },
    { selectionModeType::smCellType, "cellType" }
});


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fv::cellSetOption::setSelection(const dictionary& dict)
{
    selectionNames_.clear();

    switch (selectionMode_)
    {
        case smAll:
        {
            break;
        }
        case smGeometric:
        {
            geometricSelection_ = dict.subDict("selection");
            break;
        }
        case smPoints:
        {
            dict.readEntry("points", points_);
            break;
        }
        case smCellSet:
        {
            selectionNames_.resize(1);
            dict.readEntry("cellSet", selectionNames_.first());
            break;
        }
        case smCellZone:
        {
            if
            (
                !dict.readIfPresent("cellZones", selectionNames_)
             || selectionNames_.empty()
            )
            {
                selectionNames_.resize(1);
                dict.readEntry("cellZone", selectionNames_.first());
            }
            break;
        }
        case smCellType:
        {
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown selectionMode "
                << selectionModeTypeNames_[selectionMode_]
                << ". Valid selectionMode types : "
                << selectionModeTypeNames_
                << exit(FatalError);
        }
    }
}


void Foam::fv::cellSetOption::setVol()
{
    // Set volume information

    scalar sumVol = 0;
    for (const label celli : cells_)
    {
        sumVol += mesh_.V()[celli];
    }
    reduce(sumVol, sumOp<scalar>());

    const scalar old(V_);
    V_ = sumVol;

    // Compare volume values, stringified using current write precision
    if
    (
        Time::timeName(old, IOstream::defaultPrecision())
     != Time::timeName(V_, IOstream::defaultPrecision())
    )
    {
        Info<< indent
            << "- selected " << returnReduce(cells_.size(), sumOp<label>())
            << " cell(s) with volume " << V_ << endl;
    }
}


void Foam::fv::cellSetOption::setCellSelection()
{
    switch (selectionMode_)
    {
        case smAll:
        {
            Info<< indent << "- selecting all cells" << endl;

            cells_ = identity(mesh_.nCells());
            break;
        }
        case smGeometric:
        {
            Info<< indent << "- selecting cells geometrically" << endl;

            bitSet selectedCells
            (
                // verbosity = true
                cellBitSet::select(mesh_, geometricSelection_, true)
            );

            // From bitSet -> labels
            cells_ = selectedCells.sortedToc();
            break;
        }
        case smPoints:
        {
            Info<< indent << "- selecting cells using points" << endl;

            labelHashSet selectedCells;

            for (const point& p : points_)
            {
                const label celli = mesh_.findCell(p);

                const bool found = (celli >= 0);

                if (found)
                {
                    selectedCells.insert(celli);
                }

                if (!returnReduceOr(found))
                {
                    WarningInFunction
                        << "No owner cell found for point " << p << endl;
                }
            }

            cells_ = selectedCells.sortedToc();
            break;
        }
        case smCellSet:
        {
            Info<< indent
                << "- selecting cells using cellSet "
                << zoneName() << endl;

            cells_ = cellSet(mesh_, zoneName()).sortedToc();
            break;
        }
        case smCellZone:
        {
            Info<< indent
                << "- selecting cells using cellZones "
                << flatOutput(selectionNames_) << nl;

            const auto& zones = mesh_.cellZones();

            // Also handles groups, multiple zones etc ...
            labelList zoneIDs = zones.indices(selectionNames_);

            if (zoneIDs.empty())
            {
                FatalErrorInFunction
                    << "No matching cellZones: "
                    << flatOutput(selectionNames_) << nl
                    << "Valid zones : "
                    << flatOutput(zones.names()) << nl
                    << "Valid groups: "
                    << flatOutput(zones.groupNames())
                    << nl
                    << exit(FatalError);
            }

            if (zoneIDs.size() == 1)
            {
                cells_ = zones[zoneIDs.first()];
                // TBD: Foam::sort(cells_);
            }
            else
            {
                cells_ = zones.selection(zoneIDs).sortedToc();
            }
            break;
        }
        case smCellType:
        {
            labelHashSet selectedCells;
            const cellCellStencilObject& overlap = Stencil::New(mesh_);
            const labelList& cellTypes = overlap.cellTypes();
            forAll(cellTypes, celli)
            {
                if (cellTypes[celli] == cellCellStencil::POROUS)
                {
                    selectedCells.insert(celli);
                }
                cells_ = selectedCells.sortedToc();
            }
            break;
        }
        default:
        {
            FatalErrorInFunction
                << "Unknown selectionMode "
                << selectionModeTypeNames_[selectionMode_]
                << ". Valid selectionMode types are "
                << selectionModeTypeNames_
                << exit(FatalError);
        }
    }

    if (smAll != selectionMode_ && returnReduceAnd(cells_.empty()))
    {
        WarningInFunction
            << "No cells selected!" << endl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::cellSetOption::cellSetOption
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::option(name, modelType, dict, mesh),
    timeStart_(-1),
    duration_(0),
    selectionMode_(selectionModeTypeNames_.get("selectionMode", coeffs_)),
    selectionNames_(),
    points_(),
    geometricSelection_(),
    V_(0)
{
    Info<< incrIndent;
    read(dict);
    setSelection(coeffs_);
    setCellSelection();
    setVol();
    Info<< decrIndent;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fv::cellSetOption::isActive()
{
    if (fv::option::isActive() && inTimeLimits(mesh_.time().value()))
    {
        // Update the cell set if the mesh is changing
        if (mesh_.changing())
        {
            if (mesh_.topoChanging())
            {
                setCellSelection();
                // Force printing of new set volume
                V_ = -GREAT;
            }
            else if
            (
                selectionMode_ == smGeometric
             || selectionMode_ == smPoints
             || selectionMode_ == smCellType
            )
            {
                // Geometric selection mode(s)
                setCellSelection();
            }

            // Report new volume (if changed)
            setVol();
        }

        return true;
    }

    return false;
}


bool Foam::fv::cellSetOption::read(const dictionary& dict)
{
    if (fv::option::read(dict))
    {
        timeStart_ = -1;

        if (coeffs_.readIfPresent("timeStart", timeStart_))
        {
            coeffs_.readEntry("duration", duration_);
        }

        return true;
    }

    return false;
}


// ************************************************************************* //
