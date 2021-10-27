/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "faceSetOption.H"
#include "areaFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace fa
    {
        defineTypeNameAndDebug(faceSetOption, 0);
    }
}


const Foam::Enum
<
    Foam::fa::faceSetOption::selectionModeType
>
Foam::fa::faceSetOption::selectionModeTypeNames_
({
    { selectionModeType::smAll, "all" },
    { selectionModeType::smVolFaceZone, "volFaceZone" }
});


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fa::faceSetOption::setSelection(const dictionary& dict)
{
    switch (selectionMode_)
    {
        case smAll:
        {
            break;
        }
        case smVolFaceZone:
        {
            dict.readEntry("faceZone", faceSetName_);
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


void Foam::fa::faceSetOption::setArea()
{
    // Set area information

    scalar sumArea = 0.0;
    for (const label facei : faces_)
    {
        sumArea += regionMesh().S()[facei];
    }
    reduce(sumArea, sumOp<scalar>());

    const scalar AOld = A_;
    A_ = sumArea;

    // Convert both areas to representation using current writeprecision
    word AOldName(Time::timeName(AOld, IOstream::defaultPrecision()));
    word AName(Time::timeName(A_, IOstream::defaultPrecision()));

    if (AName != AOldName)
    {
        Info<< indent
            << "- selected " << returnReduce(faces_.size(), sumOp<label>())
            << " face(s) with area " << A_ << endl;
    }
}


void Foam::fa::faceSetOption::setFaceSelection()
{
    switch (selectionMode_)
    {
        case smVolFaceZone:
        {
            Info<< indent
                << "- selecting faces using volume-mesh faceZone "
                << faceSetName_ << endl;

            label zoneID = mesh_.faceZones().findZoneID(faceSetName_);
            if (zoneID == -1)
            {
                FatalErrorInFunction
                    << "Cannot find faceZone " << faceSetName_ << endl
                    << "Valid faceZones are " << mesh_.faceZones().names()
                    << exit(FatalError);
            }

            const faceZone& addr = mesh_.faceZones()[zoneID];

            const bitSet isZoneFace(mesh_.nFaces(), addr);

            // Do we loop over faMesh faces or over faceZone faces?
            const labelUList& faceLabels = regionMesh().faceLabels();

            label n = 0;
            for (const label facei : faceLabels)
            {
                if (isZoneFace[facei])
                {
                    n++;
                }
            }
            faces_.setSize(n);
            n = 0;
            for (const label facei : faceLabels)
            {
                if (isZoneFace[facei])
                {
                    faces_[n++] = facei;
                }
            }
            break;
        }

        case smAll:
        {
            Info<< indent << "- selecting all faces" << endl;
            faces_ = identity(regionMesh().nFaces());

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
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fa::faceSetOption::faceSetOption
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvPatch& patch
)
:
    fa::option(name, modelType, dict, patch),
    timeStart_(-1),
    duration_(0),
    selectionMode_(selectionModeTypeNames_.get("selectionMode", coeffs_)),
    faceSetName_("none"),
    A_(0)
{
    if (isActive())
    {
        Info<< incrIndent;
        read(dict);
        setSelection(coeffs_);
        setFaceSelection();
        setArea();
        Info<< decrIndent;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::fa::faceSetOption::isActive()
{
    if (fa::option::isActive() && inTimeLimits(mesh_.time().value()))
    {
        // Update the face set if the mesh is changing
        if (mesh_.changing())
        {
            if (mesh_.topoChanging())
            {
                setArea();
                // Force printing of new set area
                A_ = -GREAT;
            }

            // Report new area (if changed)
            setArea();
        }

        return true;
    }

    return false;
}


bool Foam::fa::faceSetOption::read(const dictionary& dict)
{
    if (fa::option::read(dict))
    {
        if (coeffs_.readIfPresent("timeStart", timeStart_))
        {
            coeffs_.readEntry("duration", duration_);
        }

        return true;
    }

    return false;
}


// ************************************************************************* //
