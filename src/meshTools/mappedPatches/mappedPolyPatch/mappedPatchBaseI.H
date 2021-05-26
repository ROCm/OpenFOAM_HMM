/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

inline const Foam::mappedPatchBase::sampleMode&
Foam::mappedPatchBase::mode() const
{
    return mode_;
}


inline const Foam::word& Foam::mappedPatchBase::sampleWorld() const
{
    return sampleWorld_;
}


inline const Foam::word& Foam::mappedPatchBase::sampleRegion() const
{
    if (sampleRegion_.empty())
    {
        if (!coupleGroup_.valid())
        {
            FatalErrorInFunction
                << "Supply either a regionName or a coupleGroup"
                << " for patch " << patch_.name()
                << " in region " << patch_.boundaryMesh().mesh().name()
                << exit(FatalError);
        }

        // Try and use patchGroup to find samplePatch and sampleRegion
        label samplePatchID = coupleGroup_.findOtherPatchID
        (
            patch_,
            sampleRegion_
        );

        samplePatch_ = sampleMesh().boundaryMesh()[samplePatchID].name();
    }
    return sampleRegion_;
}


inline const Foam::word& Foam::mappedPatchBase::samplePatch() const
{
    if (samplePatch_.empty())
    {
        if (!coupleGroup_.valid())
        {
            FatalErrorInFunction
                << "Supply either a patchName or a coupleGroup"
                << " for patch " << patch_.name()
                << " in region " << patch_.boundaryMesh().mesh().name()
                << exit(FatalError);
        }

        // Try and use patchGroup to find samplePatch and sampleRegion
        label samplePatchID = coupleGroup_.findOtherPatchID
        (
            patch_,
            sampleRegion_
        );

        samplePatch_ = sampleMesh().boundaryMesh()[samplePatchID].name();
    }
    return samplePatch_;
}


inline const Foam::word& Foam::mappedPatchBase::coupleGroup() const
{
    return coupleGroup_.name();
}


inline Foam::label Foam::mappedPatchBase::sampleSize() const
{
    switch (mode_)
    {
        case NEARESTPATCHFACEAMI:
        {
            return samplePolyPatch().size();
        }
        case NEARESTCELL:
        {
            return sampleMesh().nCells();
        }
        case NEARESTPATCHFACE:
        {
            return samplePolyPatch().size();
        }
        case NEARESTPATCHPOINT:
        {
            return samplePolyPatch().nPoints();
        }
        case NEARESTFACE:
        {
            const polyMesh& mesh = sampleMesh();
            return mesh.nBoundaryFaces();
        }
        default:
        {
            FatalErrorInFunction
                << "problem." << abort(FatalError);
            return -1;
        }
    }
}


inline const Foam::vector& Foam::mappedPatchBase::offset() const
{
    return offset_;
}


inline const Foam::vectorField& Foam::mappedPatchBase::offsets() const
{
    return offsets_;
}


inline Foam::label Foam::mappedPatchBase::comm() const
{
    return comm_;
}


inline bool Foam::mappedPatchBase::sameWorld() const
{
    return
        UPstream::allWorlds().size() == 1
     || UPstream::myWorld() == sampleWorld_;
}


inline bool Foam::mappedPatchBase::masterWorld() const
{
    if (sameWorld())
    {
        return true;
    }
    else
    {
        // Use ordering in allWorlds
        const label myWorld = UPstream::myWorldID();
        const label mySampleWorld =
            UPstream::allWorlds().find(sampleWorld_);
        return myWorld < mySampleWorld;
    }
}


inline bool Foam::mappedPatchBase::sameRegion() const
{
    return sameRegion_;
}


inline const Foam::mapDistribute& Foam::mappedPatchBase::map() const
{
    const polyMesh& thisMesh = patch_.boundaryMesh().mesh();
    bool topoChange =
        (sameWorld() && sampleMesh().topoChanging())
     || thisMesh.topoChanging();

    if (topoChange)
    {
        mapPtr_.reset(nullptr);
    }

    if (!mapPtr_)
    {
        calcMapping();
    }

    return *mapPtr_;
}


inline const Foam::AMIPatchToPatchInterpolation& Foam::mappedPatchBase::AMI
(
    bool forceUpdate
) const
{
    const polyMesh& thisMesh = patch_.boundaryMesh().mesh();
    bool topoChange =
        (sameWorld() && sampleMesh().topoChanging())
     || thisMesh.topoChanging();

    if (topoChange || forceUpdate)
    {
        AMIPtr_->upToDate() = false;
    }

    calcAMI();

    return *AMIPtr_;
}


// ************************************************************************* //
