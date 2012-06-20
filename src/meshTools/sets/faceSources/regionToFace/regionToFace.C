/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

#include "regionToFace.H"
#include "polyMesh.H"
#include "faceSet.H"
#include "mappedPatchBase.H"
#include "indirectPrimitivePatch.H"
#include "PatchTools.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(regionToFace, 0);

addToRunTimeSelectionTable(topoSetSource, regionToFace, word);

addToRunTimeSelectionTable(topoSetSource, regionToFace, istream);

}


Foam::topoSetSource::addToUsageTable Foam::regionToFace::usage_
(
    regionToFace::typeName,
    "\n    Usage: regionToFace <faceSet> (x y z)\n\n"
    "    Select all faces in the connected region of the faceSet"
    " starting from the point.\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Synchronise edges
void Foam::regionToFace::syncEdges
(
    const labelList& patchEdges,
    const labelList& coupledEdges,
    const bool syncNonCollocated,

    PackedBoolList& isChangedEdge,
    DynamicList<label>& changedEdges,
    labelList& allEdgeData
) const
{
    const globalMeshData& globalData = mesh_.globalData();
    const mapDistribute& map = globalData.globalEdgeSlavesMap();

    // Convert patch-edge data into cpp-edge data
    labelList cppEdgeData(map.constructSize(), labelMax);

    forAll(patchEdges, i)
    {
        label patchEdgeI = patchEdges[i];
        label coupledEdgeI = coupledEdges[i];

        if (isChangedEdge[patchEdgeI])
        {
            cppEdgeData[coupledEdgeI] = allEdgeData[patchEdgeI];
        }
    }

    // Synchronise
    globalData.syncData
    (
        cppEdgeData,
        globalData.globalEdgeSlaves(),
        (
            syncNonCollocated
          ? globalData.globalEdgeTransformedSlaves()    // transformed elems
          : labelListList(globalData.globalEdgeSlaves().size()) //no transformed
        ),
        map,
        minEqOp<label>()
    );

    // Back from cpp-edge to patch-edge data
    forAll(patchEdges, i)
    {
        label patchEdgeI = patchEdges[i];
        label coupledEdgeI = coupledEdges[i];

        if (cppEdgeData[coupledEdgeI] != labelMax)
        {
            allEdgeData[patchEdgeI] = cppEdgeData[coupledEdgeI];

            if (!isChangedEdge[patchEdgeI])
            {
                changedEdges.append(patchEdgeI);
                isChangedEdge[patchEdgeI] = true;
            }
        }
    }
}


void Foam::regionToFace::markZone
(
    const indirectPrimitivePatch& patch,
    const label procI,
    const label faceI,
    const label zoneI,
    labelList& faceZone
) const
{
    // Calculate correspondence between patch and globalData.coupledPatch.
    labelList patchEdges;
    labelList coupledEdges;
    PackedBoolList sameEdgeOrientation;
    PatchTools::matchEdges
    (
        mesh_.globalData().coupledPatch(),
        patch,

        coupledEdges,
        patchEdges,
        sameEdgeOrientation
    );


    DynamicList<label> changedEdges(patch.nEdges());
    labelList allEdgeData(patch.nEdges(), -1);
    PackedBoolList isChangedEdge(patch.nEdges());


    // Fill initial seed
    // ~~~~~~~~~~~~~~~~~

    if (Pstream::myProcNo() == procI)
    {
        const labelList& fEdges = patch.faceEdges()[faceI];
        forAll(fEdges, i)
        {
            label edgeI = fEdges[i];
            if (!isChangedEdge[edgeI])
            {
                allEdgeData[edgeI] = zoneI;
                changedEdges.append(edgeI);
                isChangedEdge[edgeI] = true;
            }
        }
    }

    syncEdges
    (
        patchEdges,
        coupledEdges,
        true,               //syncNonCollocated,

        isChangedEdge,
        changedEdges,
        allEdgeData
    );


    // Edge-Face-Edge walk across patch
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    while (true)
    {
        // From edge to face
        // ~~~~~~~~~~~~~~~~~

        DynamicList<label> changedFaces(patch.size());

        forAll(changedEdges, changedI)
        {
            label edgeI = changedEdges[changedI];
            const labelList& eFaces = patch.edgeFaces()[edgeI];

            forAll(eFaces, i)
            {
                label faceI = eFaces[i];

                if (faceZone[faceI] == -1)
                {
                    faceZone[faceI] = zoneI;
                    changedFaces.append(faceI);
                }
            }
        }


        label nChangedFaces = returnReduce(changedFaces.size(), sumOp<label>());
        if (nChangedFaces == 0)
        {
            break;
        }


        // From face to edge
        // ~~~~~~~~~~~~~~~~~

        isChangedEdge = false;
        changedEdges.clear();

        forAll(changedFaces, i)
        {
            label faceI = changedFaces[i];
            const labelList& fEdges = patch.faceEdges()[faceI];

            forAll(fEdges, fp)
            {
                label edgeI = fEdges[fp];

                if (!isChangedEdge[edgeI])
                {
                    allEdgeData[edgeI] = zoneI;
                    changedEdges.append(edgeI);
                    isChangedEdge[edgeI] = true;
                }
            }
        }

        syncEdges
        (
            patchEdges,
            coupledEdges,
            true,               //syncNonCollocated,

            isChangedEdge,
            changedEdges,
            allEdgeData
        );

        label nChangedEdges = returnReduce(changedEdges.size(), sumOp<label>());
        if (nChangedEdges == 0)
        {
            break;
        }
    }
}


void Foam::regionToFace::combine(topoSet& set, const bool add) const
{
    Info<< "    Loading subset " << setName_ << " to delimit search region."
        << endl;
    faceSet subSet(mesh_, setName_);

    indirectPrimitivePatch patch
    (
        IndirectList<face>(mesh_.faces(), subSet.toc()),
        mesh_.points()
    );

    mappedPatchBase::nearInfo ni
    (
        pointIndexHit(false, vector::zero, -1),
        Tuple2<scalar, label>
        (
            sqr(GREAT),
            Pstream::myProcNo()
        )
    );

    forAll(patch, i)
    {
        const point& fc = patch.faceCentres()[i];
        scalar d2 = magSqr(fc-nearPoint_);

        if (!ni.first().hit() || d2 < ni.second().first())
        {
            ni.second().first() = d2;
            ni.first().setHit();
            ni.first().setPoint(fc);
            ni.first().setIndex(i);
        }
    }

    // Globally reduce
    combineReduce(ni, mappedPatchBase::nearestEqOp());

    Info<< "    Found nearest face at " << ni.first().rawPoint()
        << " on processor " << ni.second().second()
        << " face " << ni.first().index()
        << " distance " << Foam::sqrt(ni.second().first()) << endl;

    labelList faceRegion(patch.size(), -1);
    markZone
    (
        patch,
        ni.second().second(),   // procI
        ni.first().index(),     // start face
        0,                      // currentZone
        faceRegion
    );

    forAll(faceRegion, faceI)
    {
        if (faceRegion[faceI] == 0)
        {
            addOrDelete(set, patch.addressing()[faceI], add);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::regionToFace::regionToFace
(
    const polyMesh& mesh,
    const word& setName,
    const point& nearPoint
)
:
    topoSetSource(mesh),
    setName_(setName),
    nearPoint_(nearPoint)
{}


// Construct from dictionary
Foam::regionToFace::regionToFace
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh),
    setName_(dict.lookup("set")),
    nearPoint_(dict.lookup("nearPoint"))
{}


// Construct from Istream
Foam::regionToFace::regionToFace
(
    const polyMesh& mesh,
    Istream& is
)
:
    topoSetSource(mesh),
    setName_(checkIs(is)),
    nearPoint_(checkIs(is))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::regionToFace::~regionToFace()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::regionToFace::applyToSet
(
    const topoSetSource::setAction action,
    topoSet& set
) const
{
    if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
    {
        Info<< "    Adding all faces of connected region of set "
            << setName_
            << " starting from point "
            << nearPoint_ << " ..." << endl;

        combine(set, true);
    }
    else if (action == topoSetSource::DELETE)
    {
        Info<< "    Removing all cells of connected region of set "
            << setName_
            << " starting from point "
            << nearPoint_ << " ..." << endl;

        combine(set, false);
    }
}


// ************************************************************************* //
