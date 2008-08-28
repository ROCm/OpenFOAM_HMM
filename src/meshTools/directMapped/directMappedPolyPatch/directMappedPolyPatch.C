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
    labelList& sampleCellProcs,
    labelList& sampleCells,
    pointField& sampleCc
) const
{
    sampleCellProcs.setSize(samples.size());
    sampleCells.setSize(samples.size());
    sampleCc.setSize(samples.size());
    sampleCc = point(-GREAT, -GREAT, -GREAT);

    {
        // Octree based search engine
        meshSearch meshSearchEngine(boundaryMesh().mesh(), false);

        forAll(samples, sampleI)
        {
            sampleCells[sampleI] = meshSearchEngine.findCell(samples[sampleI]);

            if (sampleCells[sampleI] == -1)
            {
                sampleCellProcs[sampleI] = -1;
            }
            else
            {
                sampleCellProcs[sampleI] = Pstream::myProcNo();
                sampleCc[sampleI] =
                    boundaryMesh().mesh().cellCentres()[sampleCells[sampleI]];
            }
        }
    }

    // Use only that processor that contains the sample
    Pstream::listCombineGather(sampleCells, maxEqOp<label>());
    Pstream::listCombineScatter(sampleCells);

    labelList minSampleCellProcs(sampleCellProcs);
    Pstream::listCombineGather(sampleCellProcs, maxEqOp<label>());
    Pstream::listCombineScatter(sampleCellProcs);

    Pstream::listCombineGather(sampleCc, maxEqOp<point>());
    Pstream::listCombineScatter(sampleCc);

    if (debug)
    {
        Info<< "directMappedPolyPatch::findSamples : " << endl;
        forAll(sampleCells, sampleI)
        {
            Info<< "    " << sampleI << " coord:" << samples[sampleI]
                << " found on processor:" << sampleCellProcs[sampleI]
                << " in cell:" << sampleCells[sampleI]
                << " with cc:" << sampleCc[sampleI] << endl;
        }
    }

    // Check for samples not being found
    forAll(sampleCells, sampleI)
    {
        if (sampleCells[sampleI] == -1)
        {
            FatalErrorIn
            (
                "findSamples(const pointField&, labelList&, labelList&)"
            )   << "Did not find sample " << samples[sampleI]
                << " on any processor." << exit(FatalError);
        }
    }


    // Check for samples being found in multiple cells
    {
        forAll(minSampleCellProcs, sampleI)
        {
            if (minSampleCellProcs[sampleI] == -1)
            {
                minSampleCellProcs[sampleI] = labelMax;
            }
        }
        Pstream::listCombineGather(minSampleCellProcs, minEqOp<label>());
        Pstream::listCombineScatter(minSampleCellProcs);

        forAll(minSampleCellProcs, sampleI)
        {
            if (minSampleCellProcs[sampleI] != sampleCellProcs[sampleI])
            {
                FatalErrorIn
                (
                    "directMappedPolyPatch::findSamples"
                    "(const pointField&, labelList&, labelList&)"
                )   << "Found sample " << samples[sampleI]
                    << " on processor " << minSampleCellProcs[sampleI]
                    << " and on processor " << sampleCellProcs[sampleI] << endl
                    << "This is illegal. {lease move your sampling plane."
                    << exit(FatalError);
            }
        }
    }
}

void Foam::directMappedPolyPatch::calcMapping() const
{
    if (sendCellLabelsPtr_.valid())
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

    // Find processor and cell samples are in
    labelList sampleCellProcs;
    labelList sampleCells;
    pointField sampleCc;
    findSamples(samples, sampleCellProcs, sampleCells, sampleCc);


    // Now we have all the data we need:
    // - where sample originates from (so destination when mapping):
    //   patchFaces, patchFaceProcs.
    // - cell sample is in (so source when mapping)
    //   sampleCells, sampleCellProcs.


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
            meshTools::writeOBJ(str, sampleCc[i]);
            vertI++;
            str << "l " << vertI-1 << ' ' << vertI << nl;
        }
    }


    // Check that actual offset vector (sampleCc - patchFc) is more or less
    // constant.
    if (Pstream::master())
    {
        const scalarField magOffset(mag(sampleCc - patchFc));
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
                    << " the sampled cell " << sampleCc[sampleI] << endl
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
    mapDistribute distMap(sampleCellProcs, patchFaceProcs);

    // Rework the schedule to cell data to send, face data to receive.
    schedulePtr_.reset(new List<labelPair>(distMap.schedule()));

    const labelListList& subMap = distMap.subMap();
    const labelListList& constructMap = distMap.constructMap();

    // Extract the particular data I need to send and receive.
    sendCellLabelsPtr_.reset(new labelListList(subMap.size()));
    labelListList& sendCellLabels = sendCellLabelsPtr_();

    forAll(subMap, procI)
    {
        sendCellLabels[procI] = IndirectList<label>(sampleCells, subMap[procI]);

        if (debug)
        {
            Pout<< "To proc:" << procI << " sending values of cells:"
                << sendCellLabels[procI] << endl;
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
    offset_(vector::zero),
    schedulePtr_(NULL),
    sendCellLabelsPtr_(NULL),
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
    offset_(dict.lookup("offset")),
    schedulePtr_(NULL),
    sendCellLabelsPtr_(NULL),
    receiveFaceLabelsPtr_(NULL)
{}


Foam::directMappedPolyPatch::directMappedPolyPatch
(
    const directMappedPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    polyPatch(pp, bm),
    offset_(pp.offset_),
    schedulePtr_(NULL),
    sendCellLabelsPtr_(NULL),
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
    offset_(pp.offset_),
    schedulePtr_(NULL),
    sendCellLabelsPtr_(NULL),
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
    sendCellLabelsPtr_.clear();
    receiveFaceLabelsPtr_.clear();
}


void Foam::directMappedPolyPatch::write(Ostream& os) const
{
    polyPatch::write(os);
    os.writeKeyword("offset") << offset_ << token::END_STATEMENT << nl;
}


// ************************************************************************* //
