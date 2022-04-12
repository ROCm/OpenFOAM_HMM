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

#include "binModel.H"
#include "fvMesh.H"
#include "cartesianCS.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(binModel, 0);
    defineRunTimeSelectionTable(binModel, dictionary);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<>
bool Foam::binModel::decomposePatchValues
(
    List<List<vector>>& data,
    const label bini,
    const vector& v,
    const vector& n
) const
{
    if (!decomposePatchValues_)
    {
        return false;
    }

    #ifdef FULLDEBUG
    if (data.size() != 3)
    {
        FatalErrorInFunction
            << "Inconsistent data list size - expect size 3"
            << abort(FatalError);
    }
    #endif

    data[1][bini] += n*(v & n);
    data[2][bini] += v - n*(v & n);

    return true;
}


void Foam::binModel::setCoordinateSystem
(
    const dictionary& dict,
    const word& e3Name,
    const word& e1Name
)
{
    coordSysPtr_.clear();

    if (dict.found(coordinateSystem::typeName_()))
    {
        coordSysPtr_ =
            coordinateSystem::New
            (
                mesh_,
                dict,
                coordinateSystem::typeName_()
            );

        Info<< "Setting co-ordinate system:" << nl
            << "    - type          : " << coordSysPtr_->name() << nl
            << "    - origin        : " << coordSysPtr_->origin() << nl
            << "    - e3            : " << coordSysPtr_->e3() << nl
            << "    - e1            : " << coordSysPtr_->e1() << endl;
    }
    else if (dict.found("CofR"))
    {
        const vector origin(dict.get<point>("CofR"));

        const vector e3
        (
            e3Name == word::null
          ? vector(0, 0, 1)
          : dict.get<vector>(e3Name)
        );

        const vector e1
        (
            e1Name == word::null
          ? vector(1, 0, 0)
          : dict.get<vector>(e1Name)
        );

        coordSysPtr_.reset(new coordSystem::cartesian(origin, e3, e1));
    }
    else
    {
        coordSysPtr_.reset(new coordSystem::cartesian(dict));
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::binModel::binModel
(
    const dictionary& dict,
    const fvMesh& mesh,
    const word& outputPrefix
)
:
    writeFile(mesh, outputPrefix),
    mesh_(mesh),
    decomposePatchValues_(false),
    cumulative_(false),
    coordSysPtr_(),
    nBin_(1),
    patchSet_(),
    fieldNames_(),
    cellZoneIDs_(),
    filePtrs_()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::binModel::read(const dictionary& dict)
{
    patchSet_ = mesh_.boundaryMesh().patchSet(dict.get<wordRes>("patches"));
    fieldNames_ = dict.get<wordHashSet>("fields").sortedToc();

    if (dict.found("cellZones"))
    {
        DynamicList<label> zoneIDs;
        DynamicList<wordRe> czUnmatched;
        for (const auto& cz : dict.get<wordRes>("cellZones"))
        {
            const labelList czi(mesh_.cellZones().indices(cz));

            if (czi.empty())
            {
                czUnmatched.append(cz);
            }
            else
            {
                zoneIDs.append(czi);
            }
        }

        if (czUnmatched.size())
        {
            WarningInFunction
                << "Unable to find zone(s): " << czUnmatched << nl
                << "Valid cellZones are : " << mesh_.cellZones().sortedNames()
                << endl;
        }

        cellZoneIDs_.transfer(zoneIDs);
    }

    decomposePatchValues_ = dict.get<bool>("decomposePatchValues");

    filePtrs_.setSize(fieldNames_.size());
    forAll(filePtrs_, i)
    {
        filePtrs_.set(i, createFile(fieldNames_[i] + "Bin"));
    }

    setCoordinateSystem(dict);

    return true;
}


void Foam::binModel::updateMesh(const mapPolyMesh& mpm)
{}


void Foam::binModel::movePoints(const polyMesh& mesh)
{}


// ************************************************************************* //
