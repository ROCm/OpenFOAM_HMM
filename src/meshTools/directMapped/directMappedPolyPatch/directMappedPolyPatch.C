/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
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

#include "directMappedPolyPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "ListListOps.H"
#include "meshSearch.H"
#include "mapDistribute.H"
#include "meshTools.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(directMappedPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, directMappedPolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, directMappedPolyPatch, dictionary);


    template<>
    const char* NamedEnum<directMappedPolyPatch::sampleMode, 3>::names[] =
    {
        "nearestCell",
        "nearestPatchFace",
        "nearestFace"
    };

    const NamedEnum<directMappedPolyPatch::sampleMode, 3>
        directMappedPolyPatch::sampleModeNames_;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::directMappedPolyPatch::collectSamples
(
    pointField& samples,
    labelList& patchFaceProcs,
    labelList& patchFaces,
    pointField& patchFc
) const
{

    // Collect all sample points and the faces they come from.
    List<pointField> globalFc(Pstream::nProcs());
    List<pointField> globalSamples(Pstream::nProcs());
    labelListList globalFaces(Pstream::nProcs());

    globalFc[Pstream::myProcNo()] = this->faceCentres();
    globalSamples[Pstream::myProcNo()] = globalFc[Pstream::myProcNo()]+offset_;
    globalFaces[Pstream::myProcNo()] = identity(size());

    // Distribute to all processors
    Pstream::gatherList(globalSamples);
    Pstream::scatterList(globalSamples);
    Pstream::gatherList(globalFaces);
    Pstream::scatterList(globalFaces);
    Pstream::gatherList(globalFc);
    Pstream::scatterList(globalFc);

    // Rework into straight list
    samples = ListListOps::combine<pointField>
    (
        globalSamples,
        accessOp<pointField>()
    );
    patchFaces = ListListOps::combine<labelList>
    (
        globalFaces,
        accessOp<labelList>()
    );
    patchFc = ListListOps::combine<pointField>
    (
        globalFc,
        accessOp<pointField>()
    );

    patchFaceProcs.setSize(patchFaces.size());
    labelList nPerProc
    (
        ListListOps::subSizes
        (
            globalFaces,
            accessOp<labelList>()
        )
    );
    label sampleI = 0;
    forAll(nPerProc, procI)
    {
        for (label i = 0; i < nPerProc[procI]; i++)
        {
            patchFaceProcs[sampleI++] = procI;
        }
    }
}


// Find the processor/cell containing the samples. Does not account
// for samples being found in two processors.
void Foam::directMappedPolyPatch::findSamples
(
    const pointField& samples,
    labelList& sampleProcs,
    labelList& sampleIndices,
    pointField& sampleLocations
) const
{
    const polyMesh& mesh = boundaryMesh().mesh();

    // All the info for nearest. Construct to miss
    List<nearInfo> nearest(samples.size());

    switch (mode_)
    {
        case NEARESTCELL:
        {
            if (samplePatch_.size() != 0 && samplePatch_ != "none")
            {
                FatalErrorIn
                (
                    "directMappedPolyPatch::findSamples(const pointField&,"
                    " labelList&, labelList&, pointField&) const"
                )   << "No need to supply a patch name when in "
                    << sampleModeNames_[mode_] << " mode." << exit(FatalError);
            }

            // Octree based search engine
            meshSearch meshSearchEngine(mesh, false);

            forAll(samples, sampleI)
            {
                const point& sample = samples[sampleI];

                label cellI = meshSearchEngine.findCell(sample);

                if (cellI == -1)
                {
                    nearest[sampleI].second().first() = Foam::sqr(GREAT);
                    nearest[sampleI].second().second() = Pstream::myProcNo();
                }
                else
                {
                    const point& cc = mesh.cellCentres()[cellI];

                    nearest[sampleI].first() = pointIndexHit
                    (
                        true,
                        cc,
                        cellI
                    );
                    nearest[sampleI].second().first() = magSqr(cc-sample);
                    nearest[sampleI].second().second() = Pstream::myProcNo();
                }
            }
            break;
        }

        case NEARESTPATCHFACE:
        {
            label patchI = boundaryMesh().findPatchID(samplePatch_);

            if (patchI == -1)
            {
                FatalErrorIn("directMappedPolyPatch::findSamples(..)")
                    << "Cannot find patch " << samplePatch_ << endl
                    << "Valid patches are " << boundaryMesh().names()
                    << exit(FatalError);
            }

            Random rndGen(123456);

            const polyPatch& pp = boundaryMesh()[patchI];

            if (pp.size() == 0)
            {
                forAll(samples, sampleI)
                {
                    nearest[sampleI].second().first() = Foam::sqr(GREAT);
                    nearest[sampleI].second().second() = Pstream::myProcNo();
                }
            }
            else
            {
                // patch faces
                const labelList patchFaces(identity(pp.size()) + pp.start());

                const treeBoundBox patchBb
                (
                    treeBoundBox(pp.points(), pp.meshPoints()).extend
                    (
                        rndGen,
                        1E-4
                    )
                );

                autoPtr<indexedOctree<treeDataFace> > boundaryTree
                (
                    new indexedOctree<treeDataFace>
                    (
                        treeDataFace  // all information needed to search faces
                        (
                            false,                      // do not cache bb
                            mesh,
                            patchFaces                  // boundary faces only
                        ),
                        patchBb,   // overall search domain
                        8,                              // maxLevel
                        10,                             // leafsize
                        3.0                             // duplicity
                    )
                );

                forAll(samples, sampleI)
                {
                    const point& sample = samples[sampleI];

                    pointIndexHit& nearInfo = nearest[sampleI].first();
                    nearInfo = boundaryTree().findNearest
                    (
                        sample,
                        magSqr(patchBb.span())
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
            }
            break;
        }

        case NEARESTFACE:
        {
            if (samplePatch_.size() != 0 && samplePatch_ != "none")
            {
                FatalErrorIn
                (
                    "directMappedPolyPatch::findSamples(const pointField&,"
                    " labelList&, labelList&, pointField&) const"
                )   << "No need to supply a patch name when in "
                    << sampleModeNames_[mode_] << " mode." << exit(FatalError);
            }

            // Octree based search engine
            meshSearch meshSearchEngine(mesh, false);

            forAll(samples, sampleI)
            {
                const point& sample = samples[sampleI];

                label faceI = meshSearchEngine.findNearestFace(sample);

                if (faceI == -1)
                {
                    nearest[sampleI].second().first() = Foam::sqr(GREAT);
                    nearest[sampleI].second().second() = Pstream::myProcNo();
                }
                else
                {
                    const point& fc = mesh.faceCentres()[faceI];

                    nearest[sampleI].first() = pointIndexHit
                    (
                        true,
                        fc,
                        faceI
                    );
                    nearest[sampleI].second().first() = magSqr(fc-sample);
                    nearest[sampleI].second().second() = Pstream::myProcNo();
                }
            }
            break;
        }

        default:
        {
            FatalErrorIn("directMappedPolyPatch::findSamples(..)")
                << "problem." << abort(FatalError);
        }
    }


    // Find nearest.
    Pstream::listCombineGather(nearest, nearestEqOp());
    Pstream::listCombineScatter(nearest);

    if (debug)
    {
        Info<< "directMappedPolyPatch::findSamples : " << endl;
        forAll(nearest, sampleI)
        {
            label procI = nearest[sampleI].second().second();
            label localI = nearest[sampleI].first().index();

            Info<< "    " << sampleI << " coord:"<< samples[sampleI]
                << " found on processor:" << procI
                << " in local cell/face:" << localI
                << " with cc:" << nearest[sampleI].first().rawPoint() << endl;
        }
    }

    // Check for samples not being found
    forAll(nearest, sampleI)
    {
        if (!nearest[sampleI].first().hit())
        {
            FatalErrorIn
            (
                "directMappedPolyPatch::findSamples"
                "(const pointField&, labelList&"
                ", labelList&, pointField&)"
            )   << "Did not find sample " << samples[sampleI]
                << " on any processor." << exit(FatalError);
        }
    }


    // Convert back into proc+local index
    sampleProcs.setSize(samples.size());
    sampleIndices.setSize(samples.size());
    sampleLocations.setSize(samples.size());

    forAll(nearest, sampleI)
    {
        sampleProcs[sampleI] = nearest[sampleI].second().second();
        sampleIndices[sampleI] = nearest[sampleI].first().index();
        sampleLocations[sampleI] = nearest[sampleI].first().hitPoint();
    }
}


void Foam::directMappedPolyPatch::calcMapping() const
{
    if (sendLabelsPtr_.valid())
    {
        FatalErrorIn("directMappedPolyPatch::calcMapping() const")
            << "Mapping already calculated" << exit(FatalError);
    }

    if (offset_ == vector::zero)
    {
        FatalErrorIn("directMappedPolyPatch::calcMapping() const")
            << "Invalid offset " << offset_ << endl
            << "Offset is the vector added to the patch face centres to"
            << " find the cell supplying the data."
            << exit(FatalError);
    }


    // Get global list of all samples and the processor and face they come from.
    pointField samples;
    labelList patchFaceProcs;
    labelList patchFaces;
    pointField patchFc;
    collectSamples(samples, patchFaceProcs, patchFaces, patchFc);

    // Find processor and cell/face samples are in and actual location.
    labelList sampleProcs;
    labelList sampleIndices;
    pointField sampleLocations;
    findSamples(samples, sampleProcs, sampleIndices, sampleLocations);


    // Now we have all the data we need:
    // - where sample originates from (so destination when mapping):
    //   patchFaces, patchFaceProcs.
    // - cell/face sample is in (so source when mapping)
    //   sampleIndices, sampleProcs.


    if (debug && Pstream::master())
    {
        OFstream str
        (
            boundaryMesh().mesh().time().path()
          / name()
          + "_directMapped.obj"
        );
        Pout<< "Dumping mapping as lines from patch faceCentres to"
            << " sampled cellCentres to file " << str.name() << endl;

        label vertI = 0;

        forAll(patchFc, i)
        {
            meshTools::writeOBJ(str, patchFc[i]);
            vertI++;
            meshTools::writeOBJ(str, sampleLocations[i]);
            vertI++;
            str << "l " << vertI-1 << ' ' << vertI << nl;
        }
    }


    // Check that actual offset vector (sampleLocations - patchFc) is more or
    // less constant.
    if (Pstream::master())
    {
        const scalarField magOffset(mag(sampleLocations - patchFc));
        const scalar avgOffset(average(magOffset));

        forAll(magOffset, sampleI)
        {
            if (mag(magOffset[sampleI]-avgOffset) > 0.001*avgOffset)
            {
                WarningIn("directMappedPolyPatch::calcMapping() const")
                    << "The actual cell centres picked up using offset "
                    << offset_ << " are not" << endl
                    << "    on a single plane."
                    << " This might give numerical problems." << endl
                    << "    At patchface " << patchFc[sampleI]
                    << " the sampled cell " << sampleLocations[sampleI] << endl
                    << "    is not on a plane " << avgOffset
                    << " offset from the patch." << endl
                    << "    You might want to shift your plane offset."
                    << " Set the debug flag to get a dump of sampled cells."
                    << endl;
                break;
            }
        }
    }


    // Determine schedule.
    mapDistribute distMap(sampleProcs, patchFaceProcs);

    // Rework the schedule to cell data to send, face data to receive.
    schedulePtr_.reset(new List<labelPair>(distMap.schedule()));

    const labelListList& subMap = distMap.subMap();
    const labelListList& constructMap = distMap.constructMap();

    // Extract the particular data I need to send and receive.
    sendLabelsPtr_.reset(new labelListList(subMap.size()));
    labelListList& sendLabels = sendLabelsPtr_();

    forAll(subMap, procI)
    {
        sendLabels[procI] = IndirectList<label>
        (
            sampleIndices,
            subMap[procI]
        );

        if (debug)
        {
            Pout<< "To proc:" << procI << " sending values of cells/faces:"
                << sendLabels[procI] << endl;
        }
    }

    receiveFaceLabelsPtr_.reset(new labelListList(constructMap.size()));
    labelListList& receiveFaceLabels = receiveFaceLabelsPtr_();

    forAll(constructMap, procI)
    {
        receiveFaceLabels[procI] =
            IndirectList<label>(patchFaces, constructMap[procI]);

        if (debug)
        {
            Pout<< "From proc:" << procI << " receiving values of patch faces:"
                << receiveFaceLabels[procI] << endl;
        }
    }
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::directMappedPolyPatch::directMappedPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm
)
:
    polyPatch(name, size, start, index, bm),
    mode_(NEARESTPATCHFACE),
    samplePatch_("none"),
    offset_(vector::zero),
    schedulePtr_(NULL),
    sendLabelsPtr_(NULL),
    receiveFaceLabelsPtr_(NULL)
{}


Foam::directMappedPolyPatch::directMappedPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm
)
:
    polyPatch(name, dict, index, bm),
    mode_(sampleModeNames_.read(dict.lookup("sampleMode"))),
    samplePatch_(dict.lookup("samplePatch")),
    offset_(dict.lookup("offset")),
    schedulePtr_(NULL),
    sendLabelsPtr_(NULL),
    receiveFaceLabelsPtr_(NULL)
{}


Foam::directMappedPolyPatch::directMappedPolyPatch
(
    const directMappedPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    polyPatch(pp, bm),
    mode_(pp.mode_),
    samplePatch_(pp.samplePatch_),
    offset_(pp.offset_),
    schedulePtr_(NULL),
    sendLabelsPtr_(NULL),
    receiveFaceLabelsPtr_(NULL)
{}


Foam::directMappedPolyPatch::directMappedPolyPatch
(
    const directMappedPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart
)
:
    polyPatch(pp, bm, index, newSize, newStart),
    mode_(pp.mode_),
    samplePatch_(pp.samplePatch_),
    offset_(pp.offset_),
    schedulePtr_(NULL),
    sendLabelsPtr_(NULL),
    receiveFaceLabelsPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::directMappedPolyPatch::~directMappedPolyPatch()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::directMappedPolyPatch::clearOut()
{
    schedulePtr_.clear();
    sendLabelsPtr_.clear();
    receiveFaceLabelsPtr_.clear();
}


void Foam::directMappedPolyPatch::write(Ostream& os) const
{
    polyPatch::write(os);
    os.writeKeyword("sampleMode") << sampleModeNames_[mode_]
        << token::END_STATEMENT << nl;
    os.writeKeyword("samplePatch") << samplePatch_
        << token::END_STATEMENT << nl;
    os.writeKeyword("offset") << offset_ << token::END_STATEMENT << nl;
}


// ************************************************************************* //
