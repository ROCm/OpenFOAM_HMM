/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "patchProbes.H"
#include "volFields.H"
#include "IOmanip.H"
// For 'nearInfo' helper class only
#include "mappedPatchBase.H"
#include "treeBoundBox.H"
#include "treeDataFace.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchProbes, 0);
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
        treeBoundBox overallBb(treeBoundBox::invertedBox);

        nFaces = 0;
        forAll(patchIDs, i)
        {
            const polyPatch& pp = bm[patchIDs[i]];
            forAll(pp, i)
            {
                bndFaces[nFaces++] = pp.start()+i;
                const face& f = pp[i];
                forAll(f, fp)
                {
                    const point& pt = pp.points()[f[fp]];
                    overallBb.min() = min(overallBb.min(), pt);
                    overallBb.max() = max(overallBb.max(), pt);
                }
            }
        }

        Random rndGen(123456);
        overallBb = overallBb.extend(rndGen, 1e-4);
        overallBb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
        overallBb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

        const indexedOctree<treeDataFace> boundaryTree
        (
            treeDataFace    // all information needed to search faces
            (
                false,                      // do not cache bb
                mesh,
                bndFaces                    // patch faces only
            ),
            overallBb,                      // overall search domain
            8,                              // maxLevel
            10,                             // leafsize
            3.0                             // duplicity
        );


        forAll(probeLocations(), probeI)
        {
            const point sample = probeLocations()[probeI];

            scalar span = boundaryTree.bb().mag();

            pointIndexHit info = boundaryTree.findNearest
            (
                sample,
                Foam::sqr(span)
            );

            if (!info.hit())
            {
                info = boundaryTree.findNearest(sample, Foam::sqr(GREAT));
            }

            label faceI = boundaryTree.shapes().faceLabels()[info.index()];

            const label patchi = bm.whichPatch(faceI);

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
                    faceI
                );

                sampleInfo.second().first() = magSqr(facePt-sample);
                sampleInfo.second().second() = Pstream::myProcNo();

                nearest[probeI]= sampleInfo;
            }
        }
    }


    // Find nearest.
    Pstream::listCombineGather(nearest, mappedPatchBase::nearestEqOp());
    Pstream::listCombineScatter(nearest);


    // Update actual probe locations
    forAll(nearest, sampleI)
    {
        operator[](sampleI) = nearest[sampleI].first().rawPoint();
    }


    if (debug)
    {
        Info<< "patchProbes::findElements" << " : " << endl;
        forAll(nearest, sampleI)
        {
            label procI = nearest[sampleI].second().second();
            label localI = nearest[sampleI].first().index();

            Info<< "    " << sampleI << " coord:"<< operator[](sampleI)
                << " found on processor:" << procI
                << " in local face:" << localI
                << " with location:" << nearest[sampleI].first().rawPoint()
                << endl;
        }
    }


    // Extract any local faces to sample
    elementList_.setSize(nearest.size());
    elementList_ = -1;
    faceList_.setSize(nearest.size());
    faceList_ = -1;

    forAll(nearest, sampleI)
    {
        if (nearest[sampleI].second().second() == Pstream::myProcNo())
        {
            // Store the face to sample
            faceList_[sampleI] = nearest[sampleI].first().index();
        }
    }
}


void Foam::patchProbes::readDict(const dictionary& dict)
{
    if (!dict.readIfPresent("patches", patchNames_))
    {
        word patchName(dict.lookup("patchName"));
        patchNames_ = wordReList(1, wordRe(patchName));
    }
    probes::readDict(dict);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchProbes::patchProbes
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles,
    const bool doFindElements
)
:
    probes(name, obr, dict, loadFromFiles, false)
{
    readDict(dict);

    if (doFindElements)
    {
        // Find the elements
        findElements(mesh_);

        // Open the probe streams
        prepare();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::patchProbes::~patchProbes()
{}


void Foam::patchProbes::write()
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
}


void Foam::patchProbes::read(const dictionary& dict)
{
    readDict(dict);

    // Find the elements
    findElements(mesh_);

    // Open the probe streams
    prepare();
}


// ************************************************************************* //
