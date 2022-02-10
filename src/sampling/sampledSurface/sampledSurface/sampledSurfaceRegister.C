/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

#include "sampledSurface.H"
#include "MeshedSurface.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::polySurface* Foam::sampledSurface::getRegistrySurface
(
    const objectRegistry& obr,
    word lookupName
) const
{
    if (lookupName.empty())
    {
        lookupName = this->name();
    }

    return obr.getObjectPtr<polySurface>(lookupName);
}


Foam::polySurface* Foam::sampledSurface::storeRegistrySurface
(
    objectRegistry& obr,
    word lookupName
) const
{
    if (lookupName.empty())
    {
        lookupName = this->name();
    }

    polySurface* surfptr = getRegistrySurface(obr, lookupName);

    if (!surfptr)
    {
        // Default construct and add to registry (owned by registry)
        surfptr = new polySurface(lookupName, obr, true);
    }

    surfptr->copySurface(*this);   // Copy in geometry (removes old fields)

    return surfptr;
}


bool Foam::sampledSurface::removeRegistrySurface
(
    objectRegistry& obr,
    word lookupName
) const
{
    polySurface* surfptr = getRegistrySurface(obr, lookupName);
    return obr.checkOut(surfptr);
}


// ************************************************************************* //
