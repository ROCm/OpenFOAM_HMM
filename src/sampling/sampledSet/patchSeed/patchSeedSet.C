/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

#include "patchSeedSet.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "treeBoundBox.H"
#include "treeDataFace.H"
#include "Time.H"
#include "meshTools.H"
#include "mappedPatchBase.H"
#include "indirectPrimitivePatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchSeedSet, 0);
    addToRunTimeSelectionTable(sampledSet, patchSeedSet, word);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::patchSeedSet::calcSamples
(
    DynamicList<point>& samplingPts,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces,
    DynamicList<label>& samplingSegments,
    DynamicList<scalar>& samplingCurveDist
)
{
    DebugInfo << "patchSeedSet : sampling on patches :" << endl;

    // Construct search tree for all patch faces.
    label sz = 0;
    for (const label patchi : patchSet_)
    {
        const polyPatch& pp = mesh().boundaryMesh()[patchi];

        sz += pp.size();

        DebugInfo << "    " << pp.name() << " size " << pp.size() << endl;
    }

    labelList patchFaces(sz);
    sz = 0;
    for (const label patchi : patchSet_)
    {
        const polyPatch& pp = mesh().boundaryMesh()[patchi];
        forAll(pp, i)
        {
            patchFaces[sz++] = pp.start()+i;
        }
    }


    if (!rndGenPtr_)
    {
        rndGenPtr_.reset(new Random(0));
    }
    Random& rndGen = *rndGenPtr_;


    if (selectedLocations_.size())
    {
        DynamicList<label> newPatchFaces(patchFaces.size());

        // Find the nearest patch face
        {
            // 1. All processors find nearest local patch face for all
            //    selectedLocations

            // All the info for nearest. Construct to miss
            List<mappedPatchBase::nearInfo> nearest(selectedLocations_.size());

            const indirectPrimitivePatch pp
            (
                IndirectList<face>(mesh().faces(), patchFaces),
                mesh().points()
            );

            treeBoundBox patchBb
            (
                treeBoundBox(pp.points(), pp.meshPoints()).extend
                (
                    rndGen,
                    1e-4
                )
            );
            patchBb.min() -= point::uniform(ROOTVSMALL);
            patchBb.max() += point::uniform(ROOTVSMALL);

            indexedOctree<treeDataFace> boundaryTree
            (
                treeDataFace    // all information needed to search faces
                (
                    false,      // do not cache bb
                    mesh(),
                    patchFaces  // boundary faces only
                ),
                patchBb,        // overall search domain
                8,              // maxLevel
                10,             // leafsize
                3.0             // duplicity
            );

            // Get some global dimension so all points are equally likely
            // to be found
            const scalar globalDistSqr
            (
                //magSqr
                //(
                //    boundBox
                //    (
                //        pp.points(),
                //        pp.meshPoints(),
                //        true
                //    ).span()
                //)
                GREAT
            );

            forAll(selectedLocations_, sampleI)
            {
                const point& sample = selectedLocations_[sampleI];

                pointIndexHit& nearInfo = nearest[sampleI].first();
                nearInfo = boundaryTree.findNearest
                (
                    sample,
                    globalDistSqr
                );

                if (!nearInfo.hit())
                {
                    nearest[sampleI].second().first() = Foam::sqr(GREAT);
                    nearest[sampleI].second().second() =
                        Pstream::myProcNo();
                }
                else
                {
                    point fc(pp[nearInfo.index()].centre(pp.points()));
                    nearInfo.setPoint(fc);
                    nearest[sampleI].second().first() = magSqr(fc-sample);
                    nearest[sampleI].second().second() =
                        Pstream::myProcNo();
                }
            }


            // 2. Reduce on master. Select nearest processor.

            // Find nearest. Combine on master.
            Pstream::listCombineGather(nearest, mappedPatchBase::nearestEqOp());
            Pstream::listCombineScatter(nearest);


            // 3. Pick up my local faces that have won

            forAll(nearest, sampleI)
            {
                if (nearest[sampleI].first().hit())
                {
                    label procI = nearest[sampleI].second().second();
                    label index = nearest[sampleI].first().index();

                    if (procI == Pstream::myProcNo())
                    {
                        newPatchFaces.append(pp.addressing()[index]);
                    }
                }
            }
        }

        if (debug)
        {
            Pout<< "Found " << newPatchFaces.size()
                << " out of " << selectedLocations_.size()
                << " on local processor" << endl;
        }

        patchFaces.transfer(newPatchFaces);
    }


    // Shuffle and truncate if in random mode
    label totalSize = returnReduce(patchFaces.size(), sumOp<label>());

    if (maxPoints_ < totalSize)
    {
        // Check what fraction of maxPoints_ I need to generate locally.
        label myMaxPoints =
            label(scalar(patchFaces.size())/totalSize*maxPoints_);

        labelList subset = identity(patchFaces.size());
        for (label iter = 0; iter < 4; ++iter)
        {
            forAll(subset, i)
            {
                label j = rndGen.position<label>(0, subset.size()-1);
                std::swap(subset[i], subset[j]);
            }
        }
        // Truncate
        subset.setSize(myMaxPoints);

        // Subset patchFaces
        patchFaces = labelUIndList(patchFaces, subset)();

        if (debug)
        {
            Pout<< "In random mode : selected " << patchFaces.size()
                << " faces out of " << patchFaces.size() << endl;
        }
    }


    // Get points on patchFaces.
    globalIndex globalSampleNumbers(patchFaces.size());

    samplingPts.setCapacity(patchFaces.size());
    samplingCells.setCapacity(patchFaces.size());
    samplingFaces.setCapacity(patchFaces.size());
    samplingSegments.setCapacity(patchFaces.size());
    samplingCurveDist.setCapacity(patchFaces.size());

    // For calculation of min-decomp tet base points
    (void)mesh().tetBasePtIs();

    forAll(patchFaces, i)
    {
        label facei = patchFaces[i];

        // Slightly shift point in since on warped face face-diagonal
        // decomposition might be outside cell for face-centre decomposition!
        pointIndexHit info = mappedPatchBase::facePoint
        (
            mesh(),
            facei,
            polyMesh::FACE_DIAG_TRIS
        );
        label celli = mesh().faceOwner()[facei];

        if (info.hit())
        {
            // Move the point into the cell
            const point& cc = mesh().cellCentres()[celli];
            samplingPts.append
            (
                info.hitPoint() + 1e-1*(cc-info.hitPoint())
            );
        }
        else
        {
            samplingPts.append(info.rawPoint());
        }
        samplingCells.append(celli);
        samplingFaces.append(facei);
        samplingSegments.append(0);
        samplingCurveDist.append(globalSampleNumbers.toGlobal(i));
    }
}


void Foam::patchSeedSet::genSamples()
{
    // Storage for sample points
    DynamicList<point> samplingPts;
    DynamicList<label> samplingCells;
    DynamicList<label> samplingFaces;
    DynamicList<label> samplingSegments;
    DynamicList<scalar> samplingCurveDist;

    calcSamples
    (
        samplingPts,
        samplingCells,
        samplingFaces,
        samplingSegments,
        samplingCurveDist
    );

    samplingPts.shrink();
    samplingCells.shrink();
    samplingFaces.shrink();
    samplingSegments.shrink();
    samplingCurveDist.shrink();

    // Move into *this
    setSamples
    (
        std::move(samplingPts),
        std::move(samplingCells),
        std::move(samplingFaces),
        std::move(samplingSegments),
        std::move(samplingCurveDist)
    );

    if (debug)
    {
        write(Info);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchSeedSet::patchSeedSet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const dictionary& dict
)
:
    sampledSet(name, mesh, searchEngine, dict),
    patchSet_
    (
        mesh.boundaryMesh().patchSet(dict.get<wordRes>("patches"))
    ),
    maxPoints_(dict.get<label>("maxPoints")),
    selectedLocations_
    (
        dict.getOrDefault<pointField>
        (
            "points",
            pointField(0)
        )
    )
{
    genSamples();
}


// ************************************************************************* //
