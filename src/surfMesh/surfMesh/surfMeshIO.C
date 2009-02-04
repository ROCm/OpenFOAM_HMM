/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "surfMesh.H"
#include "Time.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::surfMesh::setInstance(const fileName& inst)
{
    if (debug)
    {
        Info<< "void surfMesh::setInstance(const fileName& inst) : "
            << "Resetting file instance to " << inst << endl;
    }

    ioPoints_.writeOpt() = IOobject::AUTO_WRITE;
    ioPoints_.instance() = inst;

    ioFaces_.writeOpt() = IOobject::AUTO_WRITE;
    ioFaces_.instance() = inst;

    surfZones_.writeOpt() = IOobject::AUTO_WRITE;
    surfZones_.instance() = inst;
}


Foam::surfMesh::readUpdateState Foam::surfMesh::readUpdate()
{
    if (debug)
    {
        Info<< "surfMesh::readUpdateState surfMesh::readUpdate() : "
            << "Updating mesh based on saved data." << endl;
    }

    // Find the point and cell instance
    fileName pointsInst(time().findInstance(meshDir(), "points"));
    fileName facesInst(time().findInstance(meshDir(), "faces"));

    if (debug)
    {
        Info<< "Faces instance: old = " << facesInstance()
            << " new = " << facesInst << nl
            << "Points instance: old = " << pointsInstance()
            << " new = " << pointsInst << endl;
    }

    if (facesInst != facesInstance())
    {
        // Topological change
        if (debug)
        {
            Info << "Topological change" << endl;
        }

        clearOut();

        // Set instance to new instance. Note that points instance can differ
        // from from faces instance.
        setInstance(facesInst);
        ioPoints_.instance() = pointsInst;

        ioPoints_ = pointIOField
        (
            IOobject
            (
                "points",
                pointsInst,
                meshSubDir,
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        ioFaces_ = faceIOList
        (
            IOobject
            (
                "faces",
                facesInst,
                meshSubDir,
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        // synchronize the sizes
        MeshReference newMeshRef(ioFaces_, ioPoints_);
        MeshReference::operator=(newMeshRef);


        // Reset the surface zones
        surfZoneIOList newZones
        (
            IOobject
            (
                "surfZones",
                facesInst,
                meshSubDir,
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        // Check that zone types and names are unchanged
        bool zonesChanged = false;

        if (surfZones_.size() != newZones.size())
        {
            zonesChanged = true;
        }
        else
        {
            forAll (surfZones_, zoneI)
            {
                if (surfZones_[zoneI].name() != newZones[zoneI].name())
                {
                    zonesChanged = true;
                    break;
                }
            }
        }

        if (zonesChanged)
        {
            WarningIn("surfMesh::readUpdateState surfMesh::readUpdate()")
                << "Number of zones has changed.  This may have "
                << "unexpected consequences.  Proceed with care." << endl;

            surfZones_.transfer(newZones);
        }
        else
        {
            surfZones_.transfer(newZones);
        }

        if (zonesChanged)
        {
            return surfMesh::TOPO_PATCH_CHANGE;
        }
        else
        {
            return surfMesh::TOPO_CHANGE;
        }
    }
    else if (pointsInst != pointsInstance())
    {
        // Points moved
        if (debug)
        {
            Info << "Point motion" << endl;
        }

        clearGeom();

        ioPoints_.instance() = pointsInst;

        ioPoints_ = pointIOField
        (
            IOobject
            (
                "points",
                pointsInst,
                meshSubDir,
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        );

        return surfMesh::POINTS_MOVED;
    }
    else
    {
        if (debug)
        {
            Info << "No change" << endl;
        }

        return surfMesh::UNCHANGED;
    }
}


// ************************************************************************* //
