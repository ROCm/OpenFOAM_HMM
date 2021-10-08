/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "faMesh.H"
#include "faMeshBoundaryHalo.H"
#include "globalMeshData.H"
#include "indirectPrimitivePatch.H"
#include "edgeHashes.H"
#include "foamVtkLineWriter.H"
#include "foamVtkIndPatchWriter.H"
#include "foamVtkUIndPatchWriter.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Print out edges as point pairs
template<class PatchType>
static void printPatchEdges
(
    Ostream& os,
    const PatchType& p,
    const labelList& edgeIds,
    label maxOutput = 10
)
{
    label nOutput = 0;

    for (const label patchEdgei : edgeIds)
    {
        const edge e(p.meshEdge(patchEdgei));

        os  << "    "
            << p.points()[e.first()] << ' '
            << p.points()[e.second()] << nl;

        ++nOutput;
        if (maxOutput > 0 && nOutput >= maxOutput)
        {
            os  << " ... suppressing further output" << nl;
            break;
        }
    }
}

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::List<Foam::Pair<Foam::faMesh::patchTuple>>
Foam::faMesh::getBoundaryEdgeConnections() const
{
    const polyBoundaryMesh& pbm = mesh().boundaryMesh();

    const label nNonProcessor = pbm.nNonProcessor();

    const label nInternalEdges = patch().nInternalEdges();
    const label nBoundaryEdges = patch().nBoundaryEdges();

    // The output result:
    List<Pair<patchTuple>> bndEdgeConnections(nBoundaryEdges);

    // Map edges (mesh numbering) back to a boundary index
    EdgeMap<label> edgeToBoundaryIndex(2*nBoundaryEdges);

    label nBadEdges(0);
    labelHashSet badEdges(2*nBoundaryEdges);

    {
        // Local collection structure for accounting of patch pairs.
        // Based on 'edge' which has some hash-like insertion properties
        // that are useful.
        struct patchPairingType : public Foam::edge
        {
            label patchEdgei_ = -1;
            label meshFacei_ = -1;

            void clear()
            {
                Foam::edge::clear();
                patchEdgei_ = -1;
                meshFacei_ = -1;
            }
        };

        List<patchPairingType> patchPairings(nBoundaryEdges);

        DebugInFunction
            << "Determining required boundary edge connections, "
            << "resolving locally attached boundary edges." << endl;


        // Pass 1:
        // - setup lookup (edge -> bnd index)
        // - add owner patch for each boundary edge
        for (label bndEdgei = 0; bndEdgei < nBoundaryEdges; ++bndEdgei)
        {
            const label patchEdgei = (bndEdgei + nInternalEdges);

            edgeToBoundaryIndex.insert
            (
                patch().meshEdge(patchEdgei),
                bndEdgei
            );

            // The attached patch face. Should only be one!
            const labelList& edgeFaces = patch().edgeFaces()[patchEdgei];

            if (edgeFaces.size() != 1)
            {
                badEdges.insert(patchEdgei);
                continue;
            }

            const label patchFacei = edgeFaces[0];
            const label meshFacei = faceLabels_[patchFacei];
            const label bndFacei = (meshFacei - mesh().nInternalFaces());

            /// const label patchId = pbm.whichPatch(meshFacei);
            const label patchId = pbm.patchID()[bndFacei];

            // Primary bookkeeping
            {
                auto& tuple = bndEdgeConnections[bndEdgei].first();

                tuple.procNo(Pstream::myProcNo());
                tuple.faPatchi(patchId);  // Tag as finiteArea patch
                tuple.patchEdgei(patchEdgei);
                tuple.meshFacei(meshFacei);
            }

            // Temporary local bookkeeping
            {
                auto& pairing = patchPairings[bndEdgei];

                pairing.clear();            // invalidate
                pairing.insert(patchId);    // 'hash' into first location
                pairing.patchEdgei_ = patchEdgei;
                pairing.meshFacei_ = meshFacei;
            }
        }

        if ((nBadEdges = returnReduce(badEdges.size(), sumOp<label>())) > 0)
        {
            edgeList dumpEdges(patch().edges(), badEdges.sortedToc());

            vtk::lineWriter writer
            (
                patch().localPoints(),
                dumpEdges,
                fileName
                (
                    mesh().time().globalPath()
                  / ("faMesh-construct.nonManifoldEdges")
                )
            );

            writer.writeGeometry();

            // CellData
            writer.beginCellData();
            writer.writeProcIDs();

            FatalErrorInFunction
                << "Boundary edges not singly connected: "
                << nBadEdges << '/' << nBoundaryEdges << nl;

            printPatchEdges
            (
                FatalError,
                patch(),
                badEdges.sortedToc()
            );

            InfoInFunction
                << "(debug) wrote " << writer.output().name() << nl;

            FatalError << abort(FatalError);
        }
        badEdges.clear();

        // Pass 2:
        // Add in first connecting neighbour patch for the boundary edges.
        // Need to examine all possibly connecting (non-processor) neighbours,
        // but only need to check their boundary edges.

        label nMissing = patchPairings.size();

        for (label patchi = 0; patchi < nNonProcessor; ++patchi)
        {
            if (!nMissing) break;  // Early exit

            const polyPatch& pp = pbm[patchi];

            // Check boundary edges
            for
            (
                label patchEdgei = pp.nInternalEdges();
                patchEdgei < pp.nEdges();
                ++patchEdgei
            )
            {
                const label bndEdgei =
                    edgeToBoundaryIndex.lookup(pp.meshEdge(patchEdgei), -1);

                if (bndEdgei != -1)
                {
                    // Has a matching owner boundary edge

                    auto& pairing = patchPairings[bndEdgei];

                    // Add neighbour (patchId, patchEdgei, meshFacei)
                    // 'hash' into the second location
                    // which does not insert the same value twice
                    if (pairing.insert(patchi))
                    {
                        // The attached patch face. Should only be one!
                        const labelList& edgeFaces = pp.edgeFaces()[patchEdgei];

                        if (edgeFaces.size() != 1)
                        {
                            pairing.erase(patchi);
                            badEdges.insert(badEdges.size());
                            continue;
                        }

                        const label patchFacei = edgeFaces[0];
                        const label meshFacei = patchFacei + pp.start();

                        // The neighbour information
                        pairing.patchEdgei_ = patchEdgei;
                        pairing.meshFacei_ = meshFacei;

                        --nMissing;
                        if (!nMissing) break;  // Early exit
                    }
                }
            }
        }

        if ((nBadEdges = returnReduce(badEdges.size(), sumOp<label>())) > 0)
        {
            FatalErrorInFunction
                << "Had " << nBadEdges
                << " boundary edges with missing or multiple edge connections"
                << abort(FatalError);
        }

        // Combine local bookkeeping into final list
        badEdges.clear();
        for (label bndEdgei = 0; bndEdgei < nBoundaryEdges; ++bndEdgei)
        {
            const auto& pairing = patchPairings[bndEdgei];
            const label nbrPatchi = pairing.second();
            const label nbrPatchEdgei = pairing.patchEdgei_;
            const label nbrMeshFacei = pairing.meshFacei_;

            if (nbrMeshFacei >= 0)
            {
                // Add into primary bookkeeping
                auto& tuple = bndEdgeConnections[bndEdgei].second();

                tuple.procNo(Pstream::myProcNo());
                tuple.patchi(nbrPatchi);
                tuple.patchEdgei(nbrPatchEdgei);
                tuple.meshFacei(nbrMeshFacei);
            }
            else if (!Pstream::parRun())
            {
                badEdges.insert(nInternalEdges + bndEdgei);
            }
        }
    }


    // ~~~~~~
    // Serial - can return already
    // ~~~~~~
    if (!Pstream::parRun())
    {
        // Verbose report of missing edges - in serial

        nBadEdges = badEdges.size();
        if (nBadEdges)
        {
            edgeList dumpEdges(patch().edges(), badEdges.sortedToc());

            vtk::lineWriter writer
            (
                patch().localPoints(),
                dumpEdges,
                fileName
                (
                    mesh().time().globalPath()
                  / ("faMesh-construct.invalidEdges")
                )
            );

            writer.writeGeometry();

            // CellData
            writer.beginCellData();
            writer.writeProcIDs();

            FatalErrorInFunction
                << "Boundary edges with missing/invalid neighbours: "
                << nBadEdges << '/' << nBoundaryEdges << nl;

            printPatchEdges
            (
                FatalError,
                patch(),
                badEdges.sortedToc()
            );

            InfoInFunction
                << "(debug) wrote " << writer.output().name() << nl;

            FatalError << abort(FatalError);
        }

        // Globally consistent ordering
        patchTuple::sort(bndEdgeConnections);
        return bndEdgeConnections;
    }


    // ~~~~~~~~
    // Parallel
    // ~~~~~~~~

    DebugInFunction
        << "Creating global coupling data" << endl;

    const globalMeshData& globalData = mesh().globalData();
    const indirectPrimitivePatch& cpp = globalData.coupledPatch();
    const mapDistribute& map = globalData.globalEdgeSlavesMap();
    const label nCoupledEdges = cpp.nEdges();

    // Construct coupled edge usage with all data
    List<bool> coupledEdgesUsed(map.constructSize(), false);

    // Markup finiteArea boundary edges that are coupled across processors
    for (label cppEdgei = 0; cppEdgei < nCoupledEdges; ++cppEdgei)
    {
        coupledEdgesUsed[cppEdgei] =
            edgeToBoundaryIndex.found(cpp.meshEdge(cppEdgei));
    }

    DebugInFunction
        << "Starting sync of boundary edge topology" << endl;

    globalMeshData::syncData
    (
        coupledEdgesUsed,
        globalData.globalEdgeSlaves(),
        globalData.globalEdgeTransformedSlaves(),  // probably not used
        map,
        orEqOp<bool>()
    );

    if (debug)
    {
        label nAttached = 0;
        for (label cppEdgei = 0; cppEdgei < nCoupledEdges; ++cppEdgei)
        {
            if (coupledEdgesUsed[cppEdgei])
            {
                ++nAttached;
            }
        }

        InfoInFunction
            << "Approx "
            << returnReduce(nAttached, sumOp<label>())
            << " connected boundary edges (total, some duplicates)" << endl;
    }

    // Combine information

    // Normally 0 or 2 connections
    List<DynamicList<patchTuple, 2>> gatheredConnections(map.constructSize());

    // Map edges (mesh numbering) back to a coupled index in use
    EdgeMap<label> edgeToCoupledIndex(2*nCoupledEdges);

    // Pass 1
    // Look for attached boundary edges
    // - boundary edges from this side go into gathered connections
    // - boundary edges connected from the other side are noted for later

    for (label cppEdgei = 0; cppEdgei < nCoupledEdges; ++cppEdgei)
    {
        if (coupledEdgesUsed[cppEdgei])
        {
            const edge meshEdge(cpp.meshEdge(cppEdgei));

            const label bndEdgei =
                edgeToBoundaryIndex.lookup(meshEdge, -1);

            if (bndEdgei != -1)
            {
                // A boundary finiteEdge edge (known from this side)

                auto& gathered = gatheredConnections[cppEdgei];
                gathered.setCapacity(2);
                gathered.resize(1);
                auto& tuple = gathered.last();

                tuple = bndEdgeConnections[bndEdgei].first();
            }
            else
            {
                // Boundary edge connected from the other side
                // - mark for it to be added later

                edgeToCoupledIndex.insert(meshEdge, cppEdgei);
            }
        }
    }

    // Pass 2
    // - add previously noted boundary edges (connected from other side)
    //   into gathered connections

    badEdges.clear();
    for (label patchi = 0; patchi < nNonProcessor; ++patchi)
    {
        if (edgeToCoupledIndex.empty()) break;  // Early exit

        const polyPatch& pp = pbm[patchi];

        // Check boundary edges
        for
        (
            label patchEdgei = pp.nInternalEdges();
            patchEdgei < pp.nEdges();
            ++patchEdgei
        )
        {
            const edge meshEdge(pp.meshEdge(patchEdgei));

            const label cppEdgei =
                edgeToCoupledIndex.lookup(meshEdge, -1);

            if (cppEdgei != -1)
            {
                // A known connection

                // The attached patch face. Should only be one!
                const labelList& edgeFaces = pp.edgeFaces()[patchEdgei];

                if (edgeFaces.size() != 1)
                {
                    badEdges.insert(cppEdgei);
                    continue;
                }

                const label patchFacei = edgeFaces[0];
                const label meshFacei = patchFacei + pp.start();

                auto& gathered = gatheredConnections[cppEdgei];
                gathered.setCapacity(2);
                gathered.resize(1);
                auto& tuple = gathered.last();

                tuple.procNo(Pstream::myProcNo());
                tuple.patchi(patchi);
                tuple.patchEdgei(patchEdgei);
                tuple.meshFacei(meshFacei);

                // Do not consider again
                edgeToCoupledIndex.erase(meshEdge);

                if (edgeToCoupledIndex.empty()) break;  // Early exit
            }
        }
    }

    if ((nBadEdges = returnReduce(badEdges.size(), sumOp<label>())) > 0)
    {
        FatalErrorInFunction
            << "Had " << nBadEdges << " coupled boundary edges"
            << " with missing or multiple edge connections"
            << abort(FatalError);
    }

    DebugInFunction
        << "Starting sync of boundary edge information" << endl;

    globalMeshData::syncData
    (
        gatheredConnections,
        globalData.globalEdgeSlaves(),
        globalData.globalEdgeTransformedSlaves(),  // probably not used
        map,
        ListOps::appendEqOp<patchTuple>()
    );


    DebugInFunction
        << "Collating sync information" << endl;

    // Pick out gathered connections and add into primary bookkeeping
    for (label cppEdgei = 0; cppEdgei < nCoupledEdges; ++cppEdgei)
    {
        const auto& gathered = gatheredConnections[cppEdgei];

        const label bndEdgei =
            edgeToBoundaryIndex.lookup(cpp.meshEdge(cppEdgei), -1);

        if
        (
            // A boundary finiteEdge edge (known from this side)
            bndEdgei != -1

            // Gathered a one-to-one connection
         && gathered.size() == 2
        )
        {
            const auto& a = gathered[0];
            const auto& b = gathered[1];

            // Copy second side of connection
            auto& connection = bndEdgeConnections[bndEdgei];

            connection.second() = (connection.first() == b) ? a : b;
        }
    }

    // Check missing/invalid
    badEdges.clear();
    for (label bndEdgei = 0; bndEdgei < nBoundaryEdges; ++bndEdgei)
    {
        const auto& connection = bndEdgeConnections[bndEdgei];

        if (!connection.second().valid())
        {
            badEdges.insert(nInternalEdges + bndEdgei);
        }
    }


    if (debug & 8)
    {
        // Boundary edges
        {
            vtk::lineWriter writer
            (
                patch().localPoints(),
                patch().boundaryEdges(),
                fileName
                (
                    mesh().time().globalPath()
                  / ("faMesh-construct.boundaryEdges")
                )
            );

            writer.writeGeometry();

            // CellData
            writer.beginCellData();
            writer.writeProcIDs();

            // For each boundary edge - the associate neighbour patch
            labelList neighProc(nBoundaryEdges);
            labelList neighPatch(nBoundaryEdges);
            for (label bndEdgei = 0; bndEdgei < nBoundaryEdges; ++bndEdgei)
            {
                const auto& connection = bndEdgeConnections[bndEdgei];

                neighProc[bndEdgei] = connection.second().procNo();
                neighPatch[bndEdgei] = connection.second().patchi();
            }

            writer.write("neighProc", neighProc);
            writer.write("neighPatch", neighPatch);
        }

        // finiteArea
        {
            vtk::uindirectPatchWriter writer
            (
                patch(),
                fileName
                (
                    mesh().time().globalPath()
                  / ("faMesh-construct.faPatch")
                )
            );

            writer.writeGeometry();

            // CellData
            writer.beginCellData();
            writer.writeProcIDs();
        }
    }

    // Verbose report of missing edges
    if ((nBadEdges = returnReduce(badEdges.size(), sumOp<label>())) > 0)
    {
        edgeList dumpEdges(patch().edges(), badEdges.sortedToc());

        vtk::lineWriter writer
        (
            patch().localPoints(),
            dumpEdges,
            fileName
            (
                mesh().time().globalPath()
              / ("faMesh-construct.invalidEdges")
            )
        );

        writer.writeGeometry();

        // CellData
        writer.beginCellData();
        writer.writeProcIDs();

        FatalErrorInFunction
            << "Boundary edges with missing/invalid neighbours: "
            << nBadEdges << '/' << nBoundaryEdges << nl;

        printPatchEdges
        (
            FatalError,
            patch(),
            badEdges.sortedToc()
        );

        InfoInFunction
            << "(debug) wrote " << writer.output().name() << nl;

        FatalError << abort(FatalError);
    }


    // Globally consistent ordering
    patchTuple::sort(bndEdgeConnections);

    DebugInFunction
        << "Return sorted list of boundary connections" << endl;

    return bndEdgeConnections;
}


void Foam::faMesh::setBoundaryConnections
(
    const List<Pair<patchTuple>>& bndEdgeConnections
) const
{
    const label nInternalEdges = patch().nInternalEdges();
    const label nBoundaryEdges = patch().nBoundaryEdges();

    if (bndEdgeConnections.size() != nBoundaryEdges)
    {
        FatalErrorInFunction
            << "Sizing mismatch. Expected " << nBoundaryEdges
            << " boundary edge connections, but had "
            << bndEdgeConnections.size() << nl
            << abort(FatalError);
    }

    bndConnectPtr_.reset
    (
        new List<labelPair>(nBoundaryEdges, labelPair(-1,-1))
    );
    auto& bndConnect = *bndConnectPtr_;

    for (const auto& connection : bndEdgeConnections)
    {
        const auto& a = connection.first();
        const auto& b = connection.second();

        if (a.is_finiteArea() && a.is_localProc())
        {
            const label bndEdgei = (a.patchEdgei() - nInternalEdges);

            bndConnect[bndEdgei].first()  = b.procNo();
            bndConnect[bndEdgei].second() = b.meshFacei();
        }
        else if (b.is_finiteArea() && b.is_localProc())
        {
            const label bndEdgei = (b.patchEdgei() - nInternalEdges);

            bndConnect[bndEdgei].first()  = a.procNo();
            bndConnect[bndEdgei].second() = a.meshFacei();
        }
        else
        {
            FatalErrorInFunction
                << "Unexpected pairing input " << connection
                << " ... programming error" << nl
                << abort(FatalError);
        }
    }

    label nInvalid = 0;
    for (const auto& connection : bndConnect)
    {
        if (connection.first() < 0 || connection.second() < 0)
        {
            ++nInvalid;
        }
    }

    if (Pstream::parRun())
    {
        reduce(nInvalid, sumOp<label>());
    }

    if (nInvalid)
    {
        FatalErrorInFunction
            << "Did not properly match " << nInvalid
            << " boundary edges" << nl
            << abort(FatalError);
    }
}


void Foam::faMesh::calcBoundaryConnections() const
{
    setBoundaryConnections(this->getBoundaryEdgeConnections());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::faMesh::boundaryProcs() const
{
    const auto& connections = this->boundaryConnections();

    labelHashSet procsUsed(2*Pstream::nProcs());

    for (const labelPair& tuple : connections)
    {
        procsUsed.insert(tuple.first());
    }

    procsUsed.erase(-1);  // placeholder value
    procsUsed.erase(Pstream::myProcNo());

    return procsUsed.sortedToc();
}


Foam::List<Foam::labelPair> Foam::faMesh::boundaryProcSizes() const
{
    const auto& connections = this->boundaryConnections();

    Map<label> procCount(2*Pstream::nProcs());

    for (const labelPair& tuple : connections)
    {
        ++procCount(tuple.first());
    }
    procCount.erase(-1);  // placeholder value
    procCount.erase(Pstream::myProcNo());

    // Flatten as list
    List<labelPair> output(procCount.size());
    label count = 0;
    for (const label proci : procCount.sortedToc())
    {
        output[count].first() = proci;
        output[count].second() = procCount[proci];  // size
        ++count;
    }

    return output;
}


const Foam::faMeshBoundaryHalo& Foam::faMesh::boundaryHaloMap() const
{
    if (!haloMapPtr_)
    {
        haloMapPtr_.reset(new faMeshBoundaryHalo(*this));
    }

    return *haloMapPtr_;
}


void Foam::faMesh::calcHaloFaceGeometry() const
{
    if (haloFaceCentresPtr_ || haloFaceNormalsPtr_)
    {
        FatalErrorInFunction
            << "Halo centres/normals already calculated"
            << exit(FatalError);
    }

    DebugInFunction
        << "Calculating halo face centres/normals" << endl;

    const faceList& faces = mesh().faces();
    const pointField& points = mesh().points();

    const faMeshBoundaryHalo& halo = boundaryHaloMap();

    const labelList& inputFaceIds = halo.inputMeshFaces();

    haloFaceCentresPtr_.reset(new pointField);
    haloFaceNormalsPtr_.reset(new vectorField);

    auto& centres = *haloFaceCentresPtr_;
    auto& normals = *haloFaceNormalsPtr_;

    centres.resize(inputFaceIds.size());
    normals.resize(inputFaceIds.size());

    // My values
    forAll(inputFaceIds, i)
    {
        const face& f = faces[inputFaceIds[i]];

        centres[i] = f.centre(points);
        normals[i] = f.unitNormal(points);
    }

    // Swap information and resize
    halo.distributeSparse(centres);
    halo.distributeSparse(normals);
}


const Foam::pointField& Foam::faMesh::haloFaceCentres() const
{
    if (!haloFaceCentresPtr_ || !haloFaceNormalsPtr_)
    {
        calcHaloFaceGeometry();
    }

    return *haloFaceCentresPtr_;
}


const Foam::vectorField& Foam::faMesh::haloFaceNormals() const
{
    if (!haloFaceCentresPtr_ || !haloFaceNormalsPtr_)
    {
        calcHaloFaceGeometry();
    }

    return *haloFaceNormalsPtr_;
}


Foam::tmp<Foam::pointField>
Foam::faMesh::haloFaceCentres(const label patchi) const
{
    if (patchi < 0 || patchi >= boundary().size())
    {
        FatalErrorInFunction
            << "Patch " << patchi << " is out-of-range 0.."
            << (boundary().size()-1) << nl
            << exit(FatalError);
    }

    return tmp<pointField>::New
    (
        List<point>
        (
            boundarySubset
            (
                this->haloFaceCentres(),
                boundary()[patchi].edgeLabels()
            )
        )
    );
}


Foam::tmp<Foam::vectorField>
Foam::faMesh::haloFaceNormals(const label patchi) const
{
    if (patchi < 0 || patchi >= boundary().size())
    {
        FatalErrorInFunction
            << "Patch " << patchi << " is out-of-range 0.."
            << (boundary().size()-1) << nl
            << exit(FatalError);
    }

    return tmp<vectorField>::New
    (
        List<vector>
        (
            boundarySubset
            (
                this->haloFaceNormals(),
                boundary()[patchi].edgeLabels()
            )
        )
    );
}


// ************************************************************************* //
