/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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

#include "patchProbes.H"
#include "volFields.H"
#include "IOmanip.H"
#include "mappedPatchBase.H"
#include "treeBoundBox.H"
#include "treeDataFace.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchProbes, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        patchProbes,
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::patchProbes::findElements(const fvMesh& mesh)
{
    (void)mesh.tetBasePtIs();

    const polyBoundaryMesh& bm = mesh.boundaryMesh();

    // All the info for nearest. Construct to miss
    List<mappedPatchBase::nearInfo> nearest(this->size());

    const labelList patchIDs(bm.patchSet(patchNames_).sortedToc());

    label nFaces = 0;
    forAll(patchIDs, i)
    {
        nFaces += bm[patchIDs[i]].size();
    }

    if (nFaces > 0)
    {
        // Collect mesh faces and bounding box
        labelList bndFaces(nFaces);
        treeBoundBox overallBb(boundBox::invertedBox);

        nFaces = 0;
        forAll(patchIDs, i)
        {
            const polyPatch& pp = bm[patchIDs[i]];
            forAll(pp, i)
            {
                bndFaces[nFaces++] = pp.start()+i;
                const face& f = pp[i];

                // Without reduction.
                overallBb.add(pp.points(), f);
            }
        }

        Random rndGen(123456);
        overallBb = overallBb.extend(rndGen, 1e-4);
        overallBb.min() -= point::uniform(ROOTVSMALL);
        overallBb.max() += point::uniform(ROOTVSMALL);


        const indexedOctree<treeDataFace> boundaryTree
        (
            treeDataFace    // all information needed to search faces
            (
                false,      // do not cache bb
                mesh,
                bndFaces    // patch faces only
            ),
            overallBb,      // overall search domain
            8,              // maxLevel
            10,             // leafsize
            3.0             // duplicity
        );

        forAll(probeLocations(), probei)
        {
            const point sample = probeLocations()[probei];

            const scalar span = boundaryTree.bb().mag();

            pointIndexHit info = boundaryTree.findNearest
            (
                sample,
                Foam::sqr(span)
            );

            if (!info.hit())
            {
                info = boundaryTree.findNearest(sample, Foam::sqr(GREAT));
            }

            label facei = boundaryTree.shapes().faceLabels()[info.index()];

            const label patchi = bm.whichPatch(facei);

            if (isA<emptyPolyPatch>(bm[patchi]))
            {
                WarningInFunction
                    << " The sample point: " << sample
                    << " belongs to " << patchi
                    << " which is an empty patch. This is not permitted. "
                    << " This sample will not be included "
                    << endl;
            }
            else if (info.hit())
            {
                // Note: do we store the face centre or the actual nearest?
                // We interpolate using the faceI only though (no
                // interpolation) so it does not actually matter much, just for
                // the location written to the header.

                //const point& facePt = mesh.faceCentres()[faceI];
                const point& facePt = info.hitPoint();

                mappedPatchBase::nearInfo sampleInfo;

                sampleInfo.first() = pointIndexHit
                (
                    true,
                    facePt,
                    facei
                );

                sampleInfo.second().first() = magSqr(facePt - sample);
                sampleInfo.second().second() = Pstream::myProcNo();

                nearest[probei]= sampleInfo;
            }
        }
    }


    // Find nearest.
    Pstream::listCombineGather(nearest, mappedPatchBase::nearestEqOp());
    Pstream::listCombineScatter(nearest);

    oldPoints_.resize(this->size());

    // Update actual probe locations and store old ones
    forAll(nearest, samplei)
    {
        oldPoints_[samplei] = operator[](samplei);
        operator[](samplei) = nearest[samplei].first().rawPoint();
    }

    if (debug)
    {
        InfoInFunction << nl;
        forAll(nearest, samplei)
        {
            label proci = nearest[samplei].second().second();
            label locali = nearest[samplei].first().index();

            Info<< "    " << samplei << " coord:"<< operator[](samplei)
                << " found on processor:" << proci
                << " in local face:" << locali
                << " with location:" << nearest[samplei].first().rawPoint()
                << endl;
        }
    }

    // Extract any local faces to sample:
    // - operator[] : actual point to sample (=nearest point on patch)
    // - oldPoints_ : original provided point (might be anywhere in the mesh)
    // - elementList_   : cells, not used
    // - faceList_      : faces (now patch faces)
    // - patchIDList_   : patch corresponding to faceList
    // - processor_     : processor
    elementList_.setSize(nearest.size());
    elementList_ = -1;
    faceList_.setSize(nearest.size());
    faceList_ = -1;
    processor_.setSize(nearest.size());
    processor_ = -1;
    patchIDList_.setSize(nearest.size());
    patchIDList_ = -1;

    forAll(nearest, sampleI)
    {
        processor_[sampleI] = nearest[sampleI].second().second();

        if (nearest[sampleI].second().second() == Pstream::myProcNo())
        {
            // Store the face to sample
            faceList_[sampleI] = nearest[sampleI].first().index();
            const label facei = faceList_[sampleI];
            if (facei != -1)
            {
                processor_[sampleI] = Pstream::myProcNo();
                patchIDList_[sampleI] = bm.whichPatch(facei);
            }
        }
        reduce(processor_[sampleI], maxOp<label>());
        reduce(patchIDList_[sampleI], maxOp<label>());
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchProbes::patchProbes
(
    const word& name,
    const Time& t,
    const dictionary& dict,
    const bool loadFromFiles,
    const bool readFields
)
:
    probes(name, t, dict, loadFromFiles, false)
{
    if (readFields)
    {
        read(dict);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::patchProbes::write()
{
    if (this->size() && prepare())
    {
        sampleAndWrite(scalarFields_);
        sampleAndWrite(vectorFields_);
        sampleAndWrite(sphericalTensorFields_);
        sampleAndWrite(symmTensorFields_);
        sampleAndWrite(tensorFields_);

        sampleAndWriteSurfaceFields(surfaceScalarFields_);
        sampleAndWriteSurfaceFields(surfaceVectorFields_);
        sampleAndWriteSurfaceFields(surfaceSphericalTensorFields_);
        sampleAndWriteSurfaceFields(surfaceSymmTensorFields_);
        sampleAndWriteSurfaceFields(surfaceTensorFields_);
    }

    return true;
}


bool Foam::patchProbes::read(const dictionary& dict)
{
    if (!dict.readIfPresent("patches", patchNames_))
    {
        patchNames_.resize(1);
        patchNames_.first() = dict.get<word>("patch");
    }

    return probes::read(dict);
}


// ************************************************************************* //
