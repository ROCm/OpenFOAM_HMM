/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2018 OpenCFD Ltd.
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

#include "surfMeshSamplePlane.H"
#include "dictionary.H"
#include "polyMesh.H"
#include "volFields.H"
#include "cartesianCS.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(surfMeshSamplePlane, 0);
    addNamedToRunTimeSelectionTable
    (
        surfMeshSample,
        surfMeshSamplePlane,
        word,
        plane
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::bitSet Foam::surfMeshSamplePlane::cellSelection(const bool warn) const
{
    return cuttingPlane::cellSelection
    (
        mesh(),
        bounds_,
        zoneNames_,
        name(),
        warn
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfMeshSamplePlane::surfMeshSamplePlane
(
    const word& name,
    const polyMesh& mesh,
    const plane& planeDesc,
    const wordRes& zones,
    const bool triangulate
)
:
    surfMeshSample(name, mesh),
    SurfaceSource(planeDesc),
    zoneNames_(zones),
    bounds_(),
    triangulate_(triangulate),
    needsUpdate_(true)
{
    if (debug)
    {
        if (!zoneNames_.empty())
        {
            Info<< "cellZones " << flatOutput(zoneNames_);

            if (-1 == mesh.cellZones().findIndex(zoneNames_))
            {
                Info<< " not found!";
            }
            Info<< endl;
        }
    }
}


Foam::surfMeshSamplePlane::surfMeshSamplePlane
(
    const word& name,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    surfMeshSample(name, mesh, dict),
    SurfaceSource(plane(dict)),
    zoneNames_(),
    bounds_(dict.lookupOrDefault("bounds", boundBox::invertedBox)),
    triangulate_(dict.lookupOrDefault("triangulate", true)),
    needsUpdate_(true)
{
    if (!dict.readIfPresent("zones", zoneNames_) && dict.found("zone"))
    {
        zoneNames_.resize(1);
        dict.readEntry("zone", zoneNames_.first());
    }


    // Make plane relative to the coordinateSystem (Cartesian)
    // allow lookup from global coordinate systems
    if (dict.found(coordinateSystem::typeName_()))
    {
        coordSystem::cartesian cs
        (
            coordinateSystem::New(mesh, dict, coordinateSystem::typeName_())
        );
        plane& pln = planeDesc();

        const point  orig = cs.globalPosition(pln.origin());
        const vector norm = cs.globalVector(pln.normal());

        if (debug)
        {
            Info<< "plane " << name << " :"
                << " origin:" << origin()
                << " normal:" << normal()
                << " defined within a local coordinateSystem" << endl;
        }

        // Reassign the plane
        pln = plane(orig, norm);
    }

    if (debug)
    {
        Info<< "plane " << name << " :"
            << " origin:" << origin()
            << " normal:" << normal();

        if (!bounds_.empty())
        {
            Info<< " bounds:" << bounds_;
        }

        if (!zoneNames_.empty())
        {
            Info<< " cellZones " << flatOutput(zoneNames_);

            if (-1 == mesh.cellZones().findIndex(zoneNames_))
            {
                Info<< " not found!";
            }
        }
        Info<< endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::surfMeshSamplePlane::needsUpdate() const
{
    return needsUpdate_;
}


bool Foam::surfMeshSamplePlane::expire()
{
    // Already marked as expired
    if (needsUpdate_)
    {
        return false;
    }

    needsUpdate_ = true;
    return true;
}


bool Foam::surfMeshSamplePlane::update()
{
    if (!needsUpdate_)
    {
        return false;
    }

    performCut(mesh(), triangulate_, cellSelection(true));

    if (debug)
    {
        print(Pout);
        Pout<< endl;
    }

    // Transfer content
    getOrCreateSurfMesh().transfer
    (
        static_cast<SurfaceSource&>(*this)
    );

    needsUpdate_ = false;
    return true;
}


bool Foam::surfMeshSamplePlane::sample
(
    const word& fieldName,
    const word& sampleScheme
) const
{
    return
    (
        sampleType<scalar>(fieldName, sampleScheme)
     || sampleType<vector>(fieldName, sampleScheme)
     || sampleType<sphericalTensor>(fieldName, sampleScheme)
     || sampleType<symmTensor>(fieldName, sampleScheme)
     || sampleType<tensor>(fieldName, sampleScheme)
    );
}


void Foam::surfMeshSamplePlane::print(Ostream& os) const
{
    os  << "surfMeshSamplePlane: " << name() << " :"
        << " base:" << plane::origin()
        << " normal:" << plane::normal()
        << " triangulate:" << triangulate_
        << " faces:"  << SurfaceSource::surfFaces().size()
        << " points:" << SurfaceSource::points().size();
}


// ************************************************************************* //
