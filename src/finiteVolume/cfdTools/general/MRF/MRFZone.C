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

#include "MRFZone.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatrices.H"
#include "PackedList.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::MRFZone::MRFZone(const fvMesh& mesh, Istream& is)
:
    mesh_(mesh),
    name_(is),
    dict_(is),
    cellZoneID_(mesh_.cellZones().findZoneID(name_)),
    faceZoneID_(mesh_.faceZones().findZoneID(name_)),
    outsideFaces_(0),
    patchNames_(dict_.lookup("patches")),
    origin_(dict_.lookup("origin")),
    axis_(dict_.lookup("axis")),
    omega_(dict_.lookup("omega")),
    Omega_("Omega", omega_*axis_)
{
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    axis_ = axis_/mag(axis_);
    Omega_ = omega_*axis_;

    bool cellZoneFound = (cellZoneID_ != -1);
    reduce(cellZoneFound, orOp<bool>());

    if (!cellZoneFound)
    {
        FatalErrorIn
        (
            "Foam::MRFZone::MRFZone(const fvMesh& , const dictionary&)"
        )   << "cannot find MRF cellZone " << name_
            << exit(FatalError);
    }

    bool faceZoneFound = (faceZoneID_ != -1);
    reduce(faceZoneFound, orOp<bool>());

    if (!faceZoneFound)
    {
        // Determine faces in cell zone
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        // (does not construct cells)

        const labelList& own = mesh_.faceOwner();
        const labelList& nei = mesh_.faceNeighbour();

        // Cells in zone
        PackedBoolList zoneCell(mesh_.nCells());

        if (cellZoneID_ != -1)
        {
            const labelList& cellLabels = mesh_.cellZones()[cellZoneID_];
            forAll(cellLabels, i)
            {
                zoneCell[cellLabels[i]] = 1u;
            }
        }


        // Faces in zone
        PackedBoolList zoneFacesSet(mesh_.nFaces());

        for (label faceI = 0; faceI < mesh_.nInternalFaces(); faceI++)
        {
            if
            (
                zoneCell.get(own[faceI]) == 1u
             || zoneCell.get(nei[faceI]) == 1u
            )
            {
                zoneFacesSet[faceI] = 1u;
            }
        }
        forAll(patches, patchI)
        {
            const polyPatch& pp = patches[patchI];

            if (pp.coupled())
            {
                forAll(pp, i)
                {
                    label faceI = pp.start()+i;

                    if (zoneCell.get(own[faceI]) == 1u)
                    {
                        zoneFacesSet[faceI] = 1u;
                    }
                }
            }
        }
        syncTools::syncFaceList(mesh_, zoneFacesSet, orEqOp<unsigned int>());


        // Transfer to labelList
        label n = zoneFacesSet.count();
        outsideFaces_.setSize(n);
        n = 0;
        forAll(zoneFacesSet, faceI)
        {
            if (zoneFacesSet.get(faceI) == 1u)
            {
                outsideFaces_[n++] = faceI;
            }
        }

        Info<< nl
            << "MRFZone " << name_ << " : did not find a faceZone; using "
            << returnReduce(outsideFaces_.size(), sumOp<label>())
            << " faces internal to cellZone instead." << endl;


        // Flag use of outsideFaces
        faceZoneID_ = -2;
    }


    patchLabels_.setSize(patchNames_.size());

    forAll(patchNames_, i)
    {
        patchLabels_[i] = patches.findPatchID(patchNames_[i]);

        if (patchLabels_[i] == -1)
        {
            FatalErrorIn
            (
                "Foam::MRFZone::MRFZone(const fvMesh& , const dictionary&)"
            )   << "cannot find MRF patch " << patchNames_[i]
                << exit(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::MRFZone::addCoriolis(fvVectorMatrix& UEqn) const
{
    if (cellZoneID_ == -1)
    {
        return;
    }

    const labelList& cells = mesh_.cellZones()[cellZoneID_];
    const scalarField& V = mesh_.V();
    vectorField& Usource = UEqn.source();
    const vectorField& U = UEqn.psi();
    const vector& Omega = Omega_.value();

    forAll(cells, i)
    {
        Usource[cells[i]] -= V[cells[i]]*(Omega ^ U[cells[i]]);
    }
}


void Foam::MRFZone::relativeFlux(surfaceScalarField& phi) const
{
    if (faceZoneID_ == -1)
    {
        return;
    }

    const surfaceVectorField& Cf = mesh_.Cf();
    const surfaceVectorField& Sf = mesh_.Sf();

    const vector& origin = origin_.value();
    const vector& Omega = Omega_.value();

    // Use either zone faces or outsideFaces_
    const labelList& faces =
    (
        faceZoneID_ == -2
      ? outsideFaces_
      : mesh_.faceZones()[faceZoneID_]
    );

    forAll(faces, i)
    {
        label facei = faces[i];

        if (facei < mesh_.nInternalFaces())
        {
            phi[facei] -= (Omega ^ (Cf[facei] - origin)) & Sf[facei];
        }
        else
        {
            label patchi = mesh_.boundaryMesh().whichPatch(facei);
            label patchFacei = mesh_.boundaryMesh()[patchi].whichFace(facei);

            phi.boundaryField()[patchi][patchFacei] -=
                (Omega ^ (Cf.boundaryField()[patchi][patchFacei] - origin))
              & Sf.boundaryField()[patchi][patchFacei];
        }
    }

    forAll(patchLabels_, i)
    {
        phi.boundaryField()[patchLabels_[i]] = 0.0;
    }
}


void Foam::MRFZone::correctBoundaryVelocity(volVectorField& U) const
{
    const vector& origin = origin_.value();
    const vector& Omega = Omega_.value();

    forAll(patchLabels_, i)
    {
        U.boundaryField()[patchLabels_[i]] ==
            (Omega ^ (mesh_.Cf().boundaryField()[patchLabels_[i]] - origin));
    }
}


// ************************************************************************* //
