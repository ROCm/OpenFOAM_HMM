/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2013 OpenFOAM Foundation
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

#include "meshToMeshNew.H"
#include "OFstream.H"
#include "Time.H"
#include "globalIndex.H"
#include "mergePoints.H"
#include "treeBoundBox.H"
#include "indexedOctree.H"
#include "treeDataCell.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(meshToMeshNew, 0);

    template<>
    const char* Foam::NamedEnum
    <
        Foam::meshToMeshNew::interpolationMethod,
        3
    >::names[] =
    {
        "direct",
        "mapNearest",
        "cellVolumeWeight"
    };

    const NamedEnum<meshToMeshNew::interpolationMethod, 3>
        meshToMeshNew::interpolationMethodNames_;
}

Foam::scalar Foam::meshToMeshNew::tolerance_ = 1e-6;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::meshToMeshNew::writeConnectivity
(
    const polyMesh& src,
    const polyMesh& tgt,
    const labelListList& srcToTargetAddr
) const
{
    Pout<< "Source size = " << src.nCells() << endl;
    Pout<< "Target size = " << tgt.nCells() << endl;

    word fName("addressing_" + src.name() + "_to_" + tgt.name());

    if (Pstream::parRun())
    {
        fName = fName +  "_proc" + Foam::name(Pstream::myProcNo());
    }

    OFstream os(src.time().path()/fName + ".obj");

    label vertI = 0;
    forAll(srcToTargetAddr, i)
    {
        const labelList& tgtAddress = srcToTargetAddr[i];
        forAll(tgtAddress, j)
        {
            label tgtI = tgtAddress[j];
            const vector& c0 = src.cellCentres()[i];

            const cell& c = tgt.cells()[tgtI];
            const pointField pts(c.points(tgt.faces(), tgt.points()));
            forAll(pts, j)
            {
                const point& p = pts[j];
                os  << "v " << p.x() << ' ' << p.y() << ' ' << p.z() << nl;
                vertI++;
                os  << "v " << c0.x() << ' ' << c0.y() << ' ' << c0.z()
                    << nl;
                vertI++;
                os  << "l " << vertI - 1 << ' ' << vertI << nl;
            }
        }
    }
}


Foam::labelList Foam::meshToMeshNew::maskCells
(
    const polyMesh& src,
    const polyMesh& tgt
) const
{
    boundBox intersectBb
    (
        max(src.bounds().min(), tgt.bounds().min()),
        min(src.bounds().max(), tgt.bounds().max())
    );

    intersectBb.inflate(0.01);

    const cellList& srcCells = src.cells();
    const faceList& srcFaces = src.faces();
    const pointField& srcPts = src.points();

    DynamicList<label> cells(src.size());
    forAll(srcCells, srcI)
    {
        boundBox cellBb(srcCells[srcI].points(srcFaces, srcPts), false);
        if (intersectBb.overlaps(cellBb))
        {
            cells.append(srcI);
        }
    }

    if (debug)
    {
        Pout<< "participating source mesh cells: " << cells.size() << endl;
    }

    return cells;
}


bool Foam::meshToMeshNew::findInitialSeeds
(
    const polyMesh& src,
    const polyMesh& tgt,
    const labelList& srcCellIDs,
    const boolList& mapFlag,
    const label startSeedI,
    label& srcSeedI,
    label& tgtSeedI
) const
{
    const cellList& srcCells = src.cells();
    const faceList& srcFaces = src.faces();
    const pointField& srcPts = src.points();

    for (label i = startSeedI; i < srcCellIDs.size(); i++)
    {
        label srcI = srcCellIDs[i];

        if (mapFlag[srcI])
        {
            const pointField
                pts(srcCells[srcI].points(srcFaces, srcPts).xfer());

            switch (method_)
            {
                case imDirect:
                case imCellVolumeWeight:
                {
                    forAll(pts, ptI)
                    {
                        const point& pt = pts[ptI];
                        label tgtI = tgt.cellTree().findInside(pt);

                        if (tgtI != -1 && intersect(src, tgt, srcI, tgtI))
                        {
                            srcSeedI = srcI;
                            tgtSeedI = tgtI;

                            return true;
                        }
                    }

                    break;
                }
                case imMapNearest:
                {
                    const point& pt = pts[0];
                    pointIndexHit hit = tgt.cellTree().findNearest(pt, GREAT);

                    if (hit.hit())
                    {
                        srcSeedI = srcI;
                        tgtSeedI = hit.index();

                        return true;
                    }
                    else
                    {
                        FatalErrorIn
                        (
                            "bool Foam::meshToMeshNew::findInitialSeeds"
                            "("
                                "const polyMesh&, "
                                "const polyMesh&, "
                                "const labelList&, "
                                "const boolList&, "
                                "const label, "
                                "label&, "
                                "label&"
                            ") const"
                        )
                            << "Unable to find nearest target cell"
                            << " for source cell " << srcI
                            << " with centre "
                            << srcCells[srcI].centre(srcPts, srcFaces)
                            << abort(FatalError);
                    }

                    break;
                }
                default:
                {
                    FatalErrorIn
                    (
                        "bool Foam::meshToMeshNew::findInitialSeeds"
                        "("
                            "const polyMesh&, "
                            "const polyMesh&, "
                            "const labelList&, "
                            "const boolList&, "
                            "const label, "
                            "label&, "
                            "label&"
                        ") const"
                    )
                        << "Unhandled method: "
                        << interpolationMethodNames_[method_]
                        << abort(FatalError);
                }
            }
        }
    }

    if (debug)
    {
        Pout<< "could not find starting seed" << endl;
    }

    return false;
}


void Foam::meshToMeshNew::appendNbrCells
(
    const label cellI,
    const polyMesh& mesh,
    const DynamicList<label>& visitedCells,
    DynamicList<label>& nbrCellIDs
) const
{
    const labelList& nbrCells = mesh.cellCells()[cellI];

    // filter out cells already visited from cell neighbours
    forAll(nbrCells, i)
    {
        label nbrCellI = nbrCells[i];

        if
        (
            (findIndex(visitedCells, nbrCellI) == -1)
         && (findIndex(nbrCellIDs, nbrCellI) == -1)
        )
        {
            nbrCellIDs.append(nbrCellI);
        }
    }
}


void Foam::meshToMeshNew::normaliseWeights
(
    const word& descriptor,
    const scalarField& cellVolumes,
    const labelListList& addr,
    scalarListList& wght
) const
{
    const label nCell = returnReduce(wght.size(), sumOp<label>());

    if (nCell > 0)
    {
        scalar minW = GREAT;
        scalar maxW = -GREAT;

        forAll(wght, cellI)
        {
            scalarList& w = wght[cellI];
            scalar s = sum(w);
            scalar Vc = cellVolumes[cellI];

            forAll(w, i)
            {
                w[i] /= Vc;
            }

            minW = min(minW, s/Vc);
            maxW = max(maxW, s/Vc);
        }

        Info<< "    " << descriptor << " weights min/max = "
            << returnReduce(minW, minOp<scalar>()) << ", "
            << returnReduce(maxW, maxOp<scalar>()) << endl;
    }
}


void Foam::meshToMeshNew::calcAddressing
(
    const polyMesh& src,
    const polyMesh& tgt
)
{
    srcToTgtCellAddr_.setSize(src.nCells());
    srcToTgtCellWght_.setSize(src.nCells());

    tgtToSrcCellAddr_.setSize(tgt.nCells());
    tgtToSrcCellWght_.setSize(tgt.nCells());

    if (!src.nCells() || !tgt.nCells())
    {
        if (debug)
        {
            Pout<< "mesh interpolation: cells not on processor: Source cells = "
                << src.nCells() << ", target cells = " << tgt.nCells()
                << endl;
        }
    }

    if (!src.nCells())
    {
        return;
    }
    else if (!tgt.nCells())
    {
        if (debug)
        {
            Pout<< "mesh interpolation: hhave " << src.nCells() << " source "
                << " cells but no target cells" << endl;
        }

        return;
    }

    // (potentially) participating source mesh cells
    const labelList srcCellIDs = maskCells(src, tgt);

    // list to keep track of whether src cell can be mapped
    boolList mapFlag(src.nCells(), false);
    UIndirectList<bool>(mapFlag, srcCellIDs) = true;

    // find initial point in tgt mesh
    label srcSeedI = -1;
    label tgtSeedI = -1;
    label startSeedI = 0;

    bool startWalk =
        findInitialSeeds
        (
            src,
            tgt,
            srcCellIDs,
            mapFlag,
            startSeedI,
            srcSeedI,
            tgtSeedI
        );

    if (!startWalk)
    {
        // if meshes are collocated, after inflating the source mesh bounding
        // box tgt mesh cells may be transferred, but may still not overlap
        // with the source mesh
        return;
    }


    switch (method_)
    {
        case imDirect:
        {
            calcDirect(src, tgt, srcSeedI, tgtSeedI);
            break;
        }
        case imMapNearest:
        {
            calcMapNearest
            (
                src,
                tgt,
                srcSeedI,
                tgtSeedI,
                srcCellIDs,
                mapFlag,
                startSeedI
            );
            break;
        }
        case imCellVolumeWeight:
        {
            calcCellVolumeWeight
            (
                src,
                tgt,
                srcSeedI,
                tgtSeedI,
                srcCellIDs,
                mapFlag,
                startSeedI
            );
            break;
        }
        default:
        {
            FatalErrorIn
            (
                "void Foam::meshToMeshNew::calcAddressing"
                "("
                    "const polyMesh&, "
                    "const polyMesh&"
                ")"
            )
                << "Unknown interpolation method"
                << abort(FatalError);
        }
    }


    if (debug)
    {
        writeConnectivity(src, tgt, srcToTgtCellAddr_);
    }
}


void Foam::meshToMeshNew::calculate()
{
    Info<< "Creating mesh-to-mesh addressing for " << srcRegion_.name()
        << " and " << tgtRegion_.name() << " regions using "
        << interpolationMethodNames_[method_] << endl;

    singleMeshProc_ = calcDistribution(srcRegion_, tgtRegion_);

    if (singleMeshProc_ == -1)
    {
        // create global indexing for src and tgt meshes
        globalIndex globalSrcCells(srcRegion_.nCells());
        globalIndex globalTgtCells(tgtRegion_.nCells());

        // Create processor map of overlapping cells. This map gets
        // (possibly remote) cells from the tgt mesh such that they (together)
        // cover all of the src mesh
        autoPtr<mapDistribute> mapPtr = calcProcMap(srcRegion_, tgtRegion_);
        const mapDistribute& map = mapPtr();

        pointField newTgtPoints;
        faceList newTgtFaces;
        labelList newTgtFaceOwners;
        labelList newTgtFaceNeighbours;
        labelList newTgtCellIDs;

        distributeAndMergeCells
        (
            map,
            tgtRegion_,
            globalTgtCells,
            newTgtPoints,
            newTgtFaces,
            newTgtFaceOwners,
            newTgtFaceNeighbours,
            newTgtCellIDs
        );


        // create a new target mesh
        polyMesh newTgt
        (
            IOobject
            (
                "newTgt." + Foam::name(Pstream::myProcNo()),
                tgtRegion_.time().timeName(),
                tgtRegion_.time(),
                IOobject::NO_READ
            ),
            xferMove(newTgtPoints),
            xferMove(newTgtFaces),
            xferMove(newTgtFaceOwners),
            xferMove(newTgtFaceNeighbours),
            false                                   // no parallel comms
        );

        // create some dummy patch info
        List<polyPatch*> patches(1);
        patches[0] = new polyPatch
        (
            "defaultFaces",
            newTgt.nFaces() - newTgt.nInternalFaces(),
            newTgt.nInternalFaces(),
            0,
            newTgt.boundaryMesh(),
            word::null
        );

        newTgt.addPatches(patches);

        // force calculation of tet-base points used for point-in-cell
        (void)newTgt.tetBasePtIs();

        // force construction of cell tree
//        (void)newTgt.cellTree();

        if (debug)
        {
            Pout<< "Created newTgt mesh:" << nl
                << " old cells = " << tgtRegion_.nCells()
                << ", new cells = " << newTgt.nCells() << nl
                << " old faces = " << tgtRegion_.nFaces()
                << ", new faces = " << newTgt.nFaces() << endl;

            if (debug > 1)
            {
                Pout<< "Writing newTgt mesh: " << newTgt.name() << endl;
                newTgt.write();
            }
        }

        calcAddressing(srcRegion_, newTgt);

        // per source cell the target cell address in newTgt mesh
        forAll(srcToTgtCellAddr_, i)
        {
            labelList& addressing = srcToTgtCellAddr_[i];
            forAll(addressing, addrI)
            {
                addressing[addrI] = newTgtCellIDs[addressing[addrI]];
            }
        }

        // convert target addresses in newTgtMesh into global cell numbering
        forAll(tgtToSrcCellAddr_, i)
        {
            labelList& addressing = tgtToSrcCellAddr_[i];
            forAll(addressing, addrI)
            {
                addressing[addrI] = globalSrcCells.toGlobal(addressing[addrI]);
            }
        }

        // set up as a reverse distribute
        mapDistribute::distribute
        (
            Pstream::nonBlocking,
            List<labelPair>(),
            tgtRegion_.nCells(),
            map.constructMap(),
            map.subMap(),
            tgtToSrcCellAddr_,
            ListPlusEqOp<label>(),
            labelList()
        );

        // set up as a reverse distribute
        mapDistribute::distribute
        (
            Pstream::nonBlocking,
            List<labelPair>(),
            tgtRegion_.nCells(),
            map.constructMap(),
            map.subMap(),
            tgtToSrcCellWght_,
            ListPlusEqOp<scalar>(),
            scalarList()
        );

        // weights normalisation
        normaliseWeights
        (
            "source",
            srcRegion_.cellVolumes(),
            srcToTgtCellAddr_,
            srcToTgtCellWght_
        );

        normaliseWeights
        (
            "target",
            tgtRegion_.cellVolumes(),
            tgtToSrcCellAddr_,
            tgtToSrcCellWght_
        );

        // cache maps and reset addresses
        List<Map<label> > cMap;
        srcMapPtr_.reset
        (
            new mapDistribute(globalSrcCells, tgtToSrcCellAddr_, cMap)
        );
        tgtMapPtr_.reset
        (
            new mapDistribute(globalTgtCells, srcToTgtCellAddr_, cMap)
        );

        // collect volume intersection contributions
        reduce(V_, sumOp<scalar>());
    }
    else
    {
        calcAddressing(srcRegion_, tgtRegion_);

        normaliseWeights
        (
            "source",
            srcRegion_.cellVolumes(),
            srcToTgtCellAddr_,
            srcToTgtCellWght_
        );

        normaliseWeights
        (
            "target",
            tgtRegion_.cellVolumes(),
            tgtToSrcCellAddr_,
            tgtToSrcCellWght_
        );
    }

    Info<< "    Overlap volume: " << V_ << endl;
}


Foam::AMIPatchToPatchInterpolation::interpolationMethod
Foam::meshToMeshNew::interpolationMethodAMI
(
    const interpolationMethod method
) const
{
    switch (method_)
    {
        case imDirect:
        {
            return AMIPatchToPatchInterpolation::imDirect;
            break;
        }
        case imMapNearest:
        {
            return AMIPatchToPatchInterpolation::imMapNearest;
            break;
        }
        case imCellVolumeWeight:
        {
            return AMIPatchToPatchInterpolation::imFaceAreaWeight;
            break;
        }
        default:
        {
            FatalErrorIn
            (
                "Foam::AMIPatchToPatchInterpolation::interpolationMethod"
                "Foam::meshToMeshNew::interpolationMethodAMI"
                "("
                    "const interpolationMethod method"
                ") const"
            )
                << "Unhandled enumeration " << method_
                << abort(FatalError);
        }
    }

    return AMIPatchToPatchInterpolation::imDirect;
}


const Foam::PtrList<Foam::AMIPatchToPatchInterpolation>&
Foam::meshToMeshNew::patchAMIs() const
{
    if (patchAMIs_.empty())
    {
        patchAMIs_.setSize(srcPatchID_.size());

        forAll(srcPatchID_, i)
        {
            label srcPatchI = srcPatchID_[i];
            label tgtPatchI = tgtPatchID_[i];

            const polyPatch& srcPP = srcRegion_.boundaryMesh()[srcPatchI];
            const polyPatch& tgtPP = tgtRegion_.boundaryMesh()[tgtPatchI];

            Info<< "Creating AMI between source patch " << srcPP.name()
                << " and target patch " << tgtPP.name()
                << " using " << interpolationMethodAMI(method_)
                << endl;

            Info<< incrIndent;

            patchAMIs_.set
            (
                i,
                new AMIPatchToPatchInterpolation
                (
                    srcPP,
                    tgtPP,
                    faceAreaIntersect::tmMesh,
                    interpolationMethodAMI(method_),
                    true // flip target patch since patch normals are aligned
                )
            );

            Info<< decrIndent;
        }
    }

    return patchAMIs_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshToMeshNew::meshToMeshNew
(
    const polyMesh& src,
    const polyMesh& tgt,
    const interpolationMethod& method,
    bool interpAllPatches
)
:
    srcRegion_(src),
    tgtRegion_(tgt),
    srcPatchID_(),
    tgtPatchID_(),
    patchAMIs_(),
    srcToTgtCellAddr_(),
    tgtToSrcCellAddr_(),
    srcToTgtCellWght_(),
    tgtToSrcCellWght_(),
    method_(method),
    V_(0.0),
    singleMeshProc_(-1),
    srcMapPtr_(NULL),
    tgtMapPtr_(NULL)
{
    if (interpAllPatches)
    {
        const polyBoundaryMesh& srcBM = src.boundaryMesh();
        const polyBoundaryMesh& tgtBM = tgt.boundaryMesh();

        if (srcBM.size() != tgtBM.size())
        {
            FatalErrorIn
            (
                "Foam::meshToMeshNew::meshToMeshNew"
                "("
                    "const polyMesh&, "
                    "const polyMesh&, "
                    "const interpolationMethod&"
                ")"
            )   << "Source and target meshes are dissimiar:" << nl
                << "    Source patches: " << srcBM.size() << nl
                << "    Target patches: " << tgtBM.size() << exit(FatalError);
        }

        DynamicList<label> patchID(src.boundaryMesh().size());

        forAll(srcBM, patchI)
        {
            const polyPatch& pp = srcBM[patchI];
            if (!polyPatch::constraintType(pp.type()))
            {
                patchID.append(pp.index());
            }
        }

        srcPatchID_.transfer(patchID);
        tgtPatchID_ = srcPatchID_;
    }

    // calculate volume addressing and weights
    calculate();

    // calculate patch addressing and weights
    (void)patchAMIs();
}


Foam::meshToMeshNew::meshToMeshNew
(
    const polyMesh& src,
    const polyMesh& tgt,
    const interpolationMethod& method,
    const HashTable<word>& patchMap
)
:
    srcRegion_(src),
    tgtRegion_(tgt),
    srcPatchID_(),
    tgtPatchID_(),
    patchAMIs_(),
    srcToTgtCellAddr_(),
    tgtToSrcCellAddr_(),
    srcToTgtCellWght_(),
    tgtToSrcCellWght_(),
    method_(method),
    V_(0.0),
    singleMeshProc_(-1),
    srcMapPtr_(NULL),
    tgtMapPtr_(NULL)
{
    srcPatchID_.setSize(patchMap.size());
    tgtPatchID_.setSize(patchMap.size());

    label i = 0;
    forAllConstIter(HashTable<word>, patchMap, iter)
    {
        const word& srcPatchName = iter.key();
        const word& tgtPatchName = iter();

        const polyPatch& srcPatch = srcRegion_.boundaryMesh()[srcPatchName];
        const polyPatch& tgtPatch = tgtRegion_.boundaryMesh()[tgtPatchName];

        srcPatchID_[i] = srcPatch.index();
        tgtPatchID_[i] = tgtPatch.index();
        i++;
    }

    // calculate volume addressing and weights
    calculate();

    // calculate patch addressing and weights
    (void)patchAMIs();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::meshToMeshNew::~meshToMeshNew()
{}


// ************************************************************************* //
