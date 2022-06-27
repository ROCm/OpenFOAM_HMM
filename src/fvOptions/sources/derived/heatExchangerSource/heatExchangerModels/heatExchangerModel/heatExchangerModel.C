/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "heatExchangerModel.H"
#include "coupledPolyPatch.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(heatExchangerModel, 0);
    defineRunTimeSelectionTable(heatExchangerModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatExchangerModel::heatExchangerModel
(
    const fvMesh& mesh,
    const word& name,
    const dictionary& coeffs
)
:
    writeFile(mesh, name, typeName, coeffs),
    mesh_(mesh),
    coeffs_(coeffs),
    name_(name),
    UName_("U"),
    TName_("T"),
    phiName_("phi"),
    faceZoneName_("unknown-faceZone"),
    faceId_(),
    facePatchId_(),
    faceSign_()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::heatExchangerModel::initialise()
{
    const label zoneID = mesh_.faceZones().findZoneID(faceZoneName_);

    if (zoneID < 0)
    {
        FatalErrorInFunction
            << type() << " " << name_ << ": "
            << "    Unknown face zone name: " << faceZoneName_
            << ". Valid face zones are: " << mesh_.faceZones().names()
            << exit(FatalError);
    }

    const faceZone& fZone = mesh_.faceZones()[zoneID];

    faceId_.setSize(fZone.size());
    facePatchId_.setSize(fZone.size());
    faceSign_.setSize(fZone.size());

    label count = 0;
    forAll(fZone, i)
    {
        const label facei = fZone[i];
        label faceId = -1;
        label facePatchId = -1;
        if (mesh_.isInternalFace(facei))
        {
            faceId = facei;
            facePatchId = -1;
        }
        else
        {
            facePatchId = mesh_.boundaryMesh().whichPatch(facei);
            const polyPatch& pp = mesh_.boundaryMesh()[facePatchId];
            const auto* cpp = isA<coupledPolyPatch>(pp);

            if (cpp)
            {
                faceId = (cpp->owner() ? pp.whichFace(facei) : -1);
            }
            else if (!isA<emptyPolyPatch>(pp))
            {
                faceId = pp.whichFace(facei);
            }
            else
            {
                faceId = -1;
                facePatchId = -1;
            }
        }

        if (faceId >= 0)
        {
            if (fZone.flipMap()[i])
            {
                faceSign_[count] = -1;
            }
            else
            {
                faceSign_[count] = 1;
            }
            faceId_[count] = faceId;
            facePatchId_[count] = facePatchId;
            count++;
        }
    }
    faceId_.setSize(count);
    facePatchId_.setSize(count);
    faceSign_.setSize(count);
}


// ************************************************************************* //
