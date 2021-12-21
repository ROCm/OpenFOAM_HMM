/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
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

Application
    redistributePar

Group
    grpParallelUtilities

Description
    Redistributes existing decomposed mesh and fields according to the current
    settings in the decomposeParDict file.

    Must be run on maximum number of source and destination processors.
    Balances mesh and writes new mesh to new time directory.

    Can optionally run in decompose/reconstruct mode to decompose/reconstruct
    mesh and fields.

Usage
    \b redistributePar [OPTION]

    Options:
      - \par -decompose
        Remove any existing \a processor subdirectories and decomposes the
        mesh. Equivalent to running without processor subdirectories.

      - \par -reconstruct
        Reconstruct mesh and fields (like reconstructParMesh+reconstructPar).

      - \par -newTimes
        (in combination with -reconstruct) reconstruct only new times.

      - \par -dry-run
        (not in combination with -reconstruct) Test without actually
        decomposing.

      - \par -cellDist
        not in combination with -reconstruct) Write the cell distribution
        as a labelList, for use with 'manual'
        decomposition method and as a volScalarField for visualization.

      - \par -region \<regionName\>
        Distribute named region.

      - \par -allRegions
        Distribute all regions in regionProperties. Does not check for
        existence of processor*.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "sigFpe.H"
#include "Time.H"
#include "fvMesh.H"
#include "fvMeshTools.H"
#include "fvMeshDistribute.H"
#include "decompositionMethod.H"
#include "decompositionModel.H"
#include "timeSelector.H"
#include "PstreamReduceOps.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "IOmapDistributePolyMesh.H"
#include "IOobjectList.H"
#include "globalIndex.H"
#include "loadOrCreateMesh.H"
#include "processorFvPatchField.H"
#include "zeroGradientFvPatchFields.H"
#include "topoSet.H"
#include "regionProperties.H"
#include "basicFvGeometryScheme.H"

#include "parFvFieldReconstructor.H"
#include "parLagrangianRedistributor.H"
#include "unmappedPassivePositionParticleCloud.H"
#include "hexRef8Data.H"
#include "meshRefinement.H"
#include "pointFields.H"

#include "cyclicACMIFvPatch.H"
#include "masterUncollatedFileOperation.H"
#include "uncollatedFileOperation.H"
#include "collatedFileOperation.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const int debug(::Foam::debug::debugSwitch("redistributePar", 0));

void createTimeDirs(const fileName& path)
{
    // Get current set of local processor's time directories. Uses
    // fileHandler
    const instantList localTimeDirs(Time::findTimes(path, "constant"));

    instantList masterTimeDirs;
    if (Pstream::master())
    {
        //const bool oldParRun = Pstream::parRun(false);
        //timeDirs = Time::findTimes(path, "constant");
        //Pstream::parRun(oldParRun);  // Restore parallel state
        masterTimeDirs = localTimeDirs;
    }
    Pstream::scatter(masterTimeDirs);
    //DebugVar(masterTimeDirs);
    //DebugVar(localTimeDirs);

    // Sync any cached times (e.g. masterUncollatedFileOperation::times_)
    // since only master would have done the findTimes
    for (const instant& t : masterTimeDirs)
    {
        if (!localTimeDirs.found(t))
        {
            const fileName timePath(path/t.name());

            //Pout<< "Time:" << t << nl
            //    << "    raw       :" << timePath << nl
            //    << endl;
            mkDir(timePath);
        }
    }

    // Just to make sure remove all state and re-scan
    fileHandler().flush();
    (void)Time::findTimes(path, "constant");
}


void copyUniform
(
    const fileOperation& fh,
    const bool decompose,
    const bool reconstruct,
    const word& readTimeName,
    const objectRegistry& readDb,
    const objectRegistry& writeDb
)
{
    // Detect uniform/ at original database + time
    const IOobject readIO("uniform", readTimeName, readDb);
    const fileName readPath
    (
        fh.dirPath
        (
            false,          // local directory
            readIO,
            false           // do not search in time
        )
    );
    //if (Pstream::master() && !readPath.empty())
    if (!readPath.empty())
    {
        Info<< "Detected additional non-decomposed files in "
            << readPath << endl;

        // readPath: searching is the same for all file handlers. Typical:
        //  <case>/0.1/uniform   (parent dir, decompose mode)
        //  <case>/processor1/0.1/uniform   (redistribute/reconstruct mode)
        //  <case>/processors2/0.1/uniform  ,,
        // writePath:
        //  uncollated : <case>/0.1/uniform (reconstruct mode). Should only
        //               be done by master
        //  uncollated : <case>/processorXXX/0.1/uniform. Should be done by all.
        //  collated   : <case>/processors2/0.1/uniform. Should be done by
        //               local master only.

        // See what local directory
        const IOobject writeIO("uniform", writeDb.time().timeName(), writeDb);
        const fileName writePath
        (
            fh.objectPath
            (
                writeIO,
                word::null
            )
        );
        // Do we already have this directory?
        const fileName currentPath(fh.dirPath(false, writeIO, false));

        if (::debug)
        {
            Pout<< "    readPath   :" << readPath << endl;
            Pout<< "    writePath  :" << writePath << endl;
            Pout<< "    currentPath:" << currentPath << endl;
        }

        if (readPath == writePath)
        {
            return;
        }

        if (currentPath.empty())
        {
            if (decompose)
            {
                // All processors copy to destination
                fh.cp(readPath, writePath);
            }
            else if (reconstruct)
            {
                // Only master
                if (Pstream::master())
                {
                    const bool oldParRun = Pstream::parRun(false);
                    fh.cp(readPath, writePath);
                    Pstream::parRun(oldParRun);
                }
            }
            else
            {
                // Redistribute. If same destination path do only on master,
                // if different path do on all processors. For now check
                // if collated file handler only. tbd.
                if (isA<fileOperations::collatedFileOperation>(fh))
                {
                    // Collated
                    if (Pstream::master())
                    {
                        const bool oldParRun = Pstream::parRun(false);
                        fh.cp(readPath, writePath);
                        Pstream::parRun(oldParRun);
                    }
                }
                else
                {
                    // Assume uncollated
                    fh.cp(readPath, writePath);
                }
            }
        }
    }
}


boolList haveFacesFile(const fileName& meshPath)
{
    const fileName facesPath(meshPath/"faces");
    Info<< "Checking for mesh in " << facesPath << nl << endl;
    boolList haveMesh(Pstream::nProcs(), false);
    haveMesh[Pstream::myProcNo()] = fileHandler().isFile
    (
        fileHandler().filePath(facesPath)
    );
    Pstream::gatherList(haveMesh);
    Pstream::scatterList(haveMesh);
    Info<< "Per processor mesh availability:" << nl
        << "    " << flatOutput(haveMesh) << nl << endl;
    return haveMesh;
}


void setBasicGeometry(fvMesh& mesh)
{
    // Set the fvGeometryScheme to basic since it does not require
    // any parallel communication to construct the geometry. During
    // redistributePar one might temporarily end up with processors
    // with zero procBoundaries. Normally procBoundaries trigger geometry
    // calculation (e.g. send over cellCentres) so on the processors without
    // procBoundaries this will not happen. The call to the geometry calculation
    // is not synchronised and this might lead to a hang for geometry schemes
    // that do require synchronisation

    tmp<fvGeometryScheme> basicGeometry
    (
        fvGeometryScheme::New
        (
            mesh,
            dictionary(),
            basicFvGeometryScheme::typeName
        )
    );
    mesh.geometry(basicGeometry);
}


void printMeshData(const polyMesh& mesh)
{
    // Collect all data on master

    globalIndex globalCells(mesh.nCells());
    labelListList patchNeiProcNo(Pstream::nProcs());
    labelListList patchSize(Pstream::nProcs());
    const labelList& pPatches = mesh.globalData().processorPatches();
    patchNeiProcNo[Pstream::myProcNo()].setSize(pPatches.size());
    patchSize[Pstream::myProcNo()].setSize(pPatches.size());
    forAll(pPatches, i)
    {
        const processorPolyPatch& ppp = refCast<const processorPolyPatch>
        (
            mesh.boundaryMesh()[pPatches[i]]
        );
        patchNeiProcNo[Pstream::myProcNo()][i] = ppp.neighbProcNo();
        patchSize[Pstream::myProcNo()][i] = ppp.size();
    }
    Pstream::gatherList(patchNeiProcNo);
    Pstream::gatherList(patchSize);


    // Print stats

    globalIndex globalBoundaryFaces(mesh.nBoundaryFaces());

    label maxProcCells = 0;
    label totProcFaces = 0;
    label maxProcPatches = 0;
    label totProcPatches = 0;
    label maxProcFaces = 0;

    for (const int procI : Pstream::allProcs())
    {
        Info<< nl
            << "Processor " << procI << nl
            << "    Number of cells = " << globalCells.localSize(procI)
            << endl;

        label nProcFaces = 0;

        const labelList& nei = patchNeiProcNo[procI];

        forAll(patchNeiProcNo[procI], i)
        {
            Info<< "    Number of faces shared with processor "
                << patchNeiProcNo[procI][i] << " = " << patchSize[procI][i]
                << endl;

            nProcFaces += patchSize[procI][i];
        }

        Info<< "    Number of processor patches = " << nei.size() << nl
            << "    Number of processor faces = " << nProcFaces << nl
            << "    Number of boundary faces = "
            << globalBoundaryFaces.localSize(procI)-nProcFaces << endl;

        maxProcCells = max(maxProcCells, globalCells.localSize(procI));
        totProcFaces += nProcFaces;
        totProcPatches += nei.size();
        maxProcPatches = max(maxProcPatches, nei.size());
        maxProcFaces = max(maxProcFaces, nProcFaces);
    }

    // Stats

    scalar avgProcCells = scalar(globalCells.size())/Pstream::nProcs();
    scalar avgProcPatches = scalar(totProcPatches)/Pstream::nProcs();
    scalar avgProcFaces = scalar(totProcFaces)/Pstream::nProcs();

    // In case of all faces on one processor. Just to avoid division by 0.
    if (totProcPatches == 0)
    {
        avgProcPatches = 1;
    }
    if (totProcFaces == 0)
    {
        avgProcFaces = 1;
    }

    Info<< nl
        << "Number of processor faces = " << totProcFaces/2 << nl
        << "Max number of cells = " << maxProcCells
        << " (" << 100.0*(maxProcCells-avgProcCells)/avgProcCells
        << "% above average " << avgProcCells << ")" << nl
        << "Max number of processor patches = " << maxProcPatches
        << " (" << 100.0*(maxProcPatches-avgProcPatches)/avgProcPatches
        << "% above average " << avgProcPatches << ")" << nl
        << "Max number of faces between processors = " << maxProcFaces
        << " (" << 100.0*(maxProcFaces-avgProcFaces)/avgProcFaces
        << "% above average " << avgProcFaces << ")" << nl
        << endl;
}


// Debugging: write volScalarField with decomposition for post processing.
void writeDecomposition
(
    const word& name,
    const fvMesh& mesh,
    const labelUList& decomp
)
{
    // Write the decomposition as labelList for use with 'manual'
    // decomposition method.
    labelIOList cellDecomposition
    (
        IOobject
        (
            "cellDecomposition",
            mesh.facesInstance(),  // mesh read from facesInstance
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        decomp
    );
    cellDecomposition.write();

    Info<< "Writing wanted cell distribution to volScalarField " << name
        << " for postprocessing purposes." << nl << endl;

    volScalarField procCells
    (
        IOobject
        (
            name,
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            false                   // do not register
        ),
        mesh,
        dimensionedScalar(name, dimless, -1),
        zeroGradientFvPatchScalarField::typeName
    );

    forAll(procCells, cI)
    {
        procCells[cI] = decomp[cI];
    }

    procCells.correctBoundaryConditions();
    procCells.write();
}


void determineDecomposition
(
    const Time& baseRunTime,
    const fileName& decompDictFile, // optional location for decomposeParDict
    const bool decompose,       // decompose, i.e. read from undecomposed case
    const fileName& proc0CaseName,
    const fvMesh& mesh,
    const bool writeCellDist,

    label& nDestProcs,
    labelList& decomp
)
{
    // Read decomposeParDict (on all processors)
    const decompositionModel& method = decompositionModel::New
    (
        mesh,
        decompDictFile
    );

    decompositionMethod& decomposer = method.decomposer();

    if (!decomposer.parallelAware())
    {
        WarningInFunction
            << "You have selected decomposition method "
            << decomposer.typeName
            << " which does" << nl
            << "not synchronise the decomposition across"
            << " processor patches." << nl
            << "    You might want to select a decomposition method"
            << " which is aware of this. Continuing."
            << endl;
    }

    Time& tm = const_cast<Time&>(mesh.time());

    const bool oldProcCase = tm.processorCase();
    if (Pstream::master() && decompose)
    {
        Info<< "Setting caseName to " << baseRunTime.caseName()
            << " to read decomposeParDict" << endl;
        tm.caseName() = baseRunTime.caseName();
        tm.processorCase(false);
    }

    scalarField cellWeights;
    if (method.found("weightField"))
    {
        word weightName = method.get<word>("weightField");

        volScalarField weights
        (
            IOobject
            (
                weightName,
                tm.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );
        cellWeights = weights.internalField();
    }

    nDestProcs = decomposer.nDomains();
    decomp = decomposer.decompose(mesh, cellWeights);

    if (Pstream::master() && decompose)
    {
        Info<< "Restoring caseName to " << proc0CaseName << endl;
        tm.caseName() = proc0CaseName;
        tm.processorCase(oldProcCase);
    }

    // Dump decomposition to volScalarField
    if (writeCellDist)
    {
        // Note: on master make sure to write to processor0
        if (decompose)
        {
            if (Pstream::master())
            {
                const bool oldParRun = Pstream::parRun(false);

                Info<< "Setting caseName to " << baseRunTime.caseName()
                    << " to write undecomposed cellDist" << endl;

                tm.caseName() = baseRunTime.caseName();
                tm.processorCase(false);
                writeDecomposition("cellDist", mesh, decomp);
                Info<< "Restoring caseName to " << proc0CaseName << endl;
                tm.caseName() = proc0CaseName;
                tm.processorCase(oldProcCase);

                Pstream::parRun(oldParRun);
            }
        }
        else
        {
            writeDecomposition("cellDist", mesh, decomp);
        }
    }
}


// Write addressing if decomposing (1 to many) or reconstructing (many to 1)
void writeProcAddressing
(
    autoPtr<fileOperation>&& writeHandler,
    const fvMesh& mesh,
    const mapDistributePolyMesh& map,
    const bool decompose
)
{
    Info<< "Writing procAddressing files to " << mesh.facesInstance()
        << endl;

    labelIOList cellMap
    (
        IOobject
        (
            "cellProcAddressing",
            mesh.facesInstance(),
            polyMesh::meshSubDir,
            mesh,
            IOobject::NO_READ
        ),
        0
    );

    labelIOList faceMap
    (
        IOobject
        (
            "faceProcAddressing",
            mesh.facesInstance(),
            polyMesh::meshSubDir,
            mesh,
            IOobject::NO_READ
        ),
        0
    );

    labelIOList pointMap
    (
        IOobject
        (
            "pointProcAddressing",
            mesh.facesInstance(),
            polyMesh::meshSubDir,
            mesh,
            IOobject::NO_READ
        ),
        0
    );

    labelIOList patchMap
    (
        IOobject
        (
            "boundaryProcAddressing",
            mesh.facesInstance(),
            polyMesh::meshSubDir,
            mesh,
            IOobject::NO_READ
        ),
        0
    );

    // Decomposing: see how cells moved from undecomposed case
    if (decompose)
    {
        cellMap = identity(map.nOldCells());
        map.distributeCellData(cellMap);

        faceMap = identity(map.nOldFaces());
        {
            const mapDistribute& faceDistMap = map.faceMap();

            if (faceDistMap.subHasFlip() || faceDistMap.constructHasFlip())
            {
                // Offset by 1
                faceMap = faceMap + 1;
            }
            // Apply face flips
            mapDistributeBase::distribute
            (
                Pstream::commsTypes::nonBlocking,
                List<labelPair>(),
                faceDistMap.constructSize(),
                faceDistMap.subMap(),
                faceDistMap.subHasFlip(),
                faceDistMap.constructMap(),
                faceDistMap.constructHasFlip(),
                faceMap,
                flipLabelOp()
            );
        }

        pointMap = identity(map.nOldPoints());
        map.distributePointData(pointMap);

        patchMap = identity(map.oldPatchSizes().size());
        const mapDistribute& patchDistMap = map.patchMap();
        // Use explicit distribute since we need to provide a null value
        // (for new patches) and this is the only call that allow us to
        // provide one ...
        mapDistributeBase::distribute
        (
            Pstream::commsTypes::nonBlocking,
            List<labelPair>(),
            patchDistMap.constructSize(),
            patchDistMap.subMap(),
            patchDistMap.subHasFlip(),
            patchDistMap.constructMap(),
            patchDistMap.constructHasFlip(),
            patchMap,
            label(-1),
            eqOp<label>(),
            flipOp(),
            UPstream::msgType()
        );
    }
    else    // reconstruct
    {
        cellMap = identity(mesh.nCells());
        map.cellMap().reverseDistribute(map.nOldCells(), cellMap);

        faceMap = identity(mesh.nFaces());
        {
            const mapDistribute& faceDistMap = map.faceMap();

            if (faceDistMap.subHasFlip() || faceDistMap.constructHasFlip())
            {
                // Offset by 1
                faceMap = faceMap + 1;
            }

            mapDistributeBase::distribute
            (
                Pstream::commsTypes::nonBlocking,
                List<labelPair>(),
                map.nOldFaces(),
                faceDistMap.constructMap(),
                faceDistMap.constructHasFlip(),
                faceDistMap.subMap(),
                faceDistMap.subHasFlip(),
                faceMap,
                flipLabelOp()
            );
        }

        pointMap = identity(mesh.nPoints());
        map.pointMap().reverseDistribute(map.nOldPoints(), pointMap);

        const mapDistribute& patchDistMap = map.patchMap();
        patchMap = identity(mesh.boundaryMesh().size());
        patchDistMap.reverseDistribute
        (
            map.oldPatchSizes().size(),
            label(-1),
            patchMap
        );
    }

    autoPtr<fileOperation> defaultHandler;
    if (writeHandler.valid())
    {
        defaultHandler = fileHandler(std::move(writeHandler));
    }

    const bool cellOk = cellMap.write();
    const bool faceOk = faceMap.write();
    const bool pointOk = pointMap.write();
    const bool patchOk = patchMap.write();

    if (defaultHandler.valid())
    {
        writeHandler = fileHandler(std::move(defaultHandler));
    }

    if (!cellOk || !faceOk || !pointOk || !patchOk)
    {
        WarningInFunction
            << "Failed to write " << cellMap.objectPath()
            << ", " << faceMap.objectPath()
            << ", " << pointMap.objectPath()
            << ", " << patchMap.objectPath()
            << endl;
    }
}


// Remove addressing
void removeProcAddressing(const polyMesh& mesh)
{
    for (const auto prefix : {"boundary", "cell", "face", "point"})
    {
        IOobject io
        (
            prefix + word("ProcAddressing"),
            mesh.facesInstance(),
            polyMesh::meshSubDir,
            mesh
        );

        const fileName procFile(io.objectPath());
        rm(procFile);
    }
}


// Generic mesh-based field reading
template<class GeoField>
void readField
(
    const IOobject& io,
    const fvMesh& mesh,
    const label i,
    PtrList<GeoField>& fields
)
{
    fields.set(i, new GeoField(io, mesh));
}


// Definition of readField for GeometricFields only
template<class Type, template<class> class PatchField, class GeoMesh>
void readField
(
    const IOobject& io,
    const fvMesh& mesh,
    const label i,
    PtrList<GeometricField<Type, PatchField, GeoMesh>>& fields
)
{
    fields.set
    (
        i,
        new GeometricField<Type, PatchField, GeoMesh>(io, mesh, false)
    );
}


// Read vol or surface fields
template<class GeoField>
void readFields
(
    const boolList& haveMesh,
    const fvMesh& mesh,
    const autoPtr<fvMeshSubset>& subsetterPtr,
    IOobjectList& allObjects,
    PtrList<GeoField>& fields
)
{
    // Get my objects of type
    IOobjectList objects(allObjects.lookupClass(GeoField::typeName));

    // Check that we all have all objects
    wordList objectNames = objects.sortedNames();

    // Get master names
    wordList masterNames(objectNames);
    Pstream::scatter(masterNames);

    if (haveMesh[Pstream::myProcNo()] && objectNames != masterNames)
    {
        FatalErrorInFunction
            << "Objects not synchronised across processors." << nl
            << "Master has " << flatOutput(masterNames) << nl
            << "Processor " << Pstream::myProcNo()
            << " has " << flatOutput(objectNames)
            << exit(FatalError);
    }

    fields.setSize(masterNames.size());

    // Have master send all fields to processors that don't have a mesh
    if (Pstream::master())
    {
        forAll(masterNames, i)
        {
            const word& name = masterNames[i];
            IOobject& io = *objects[name];
            io.writeOpt(IOobject::AUTO_WRITE);

            // Load field (but not oldTime)
            //const bool oldParRun = Pstream::parRun(false);
            readField(io, mesh, i, fields);
            //Pstream::parRun(oldParRun);

            // Create zero sized field and send
            if (subsetterPtr)
            {
                const bool oldParRun = Pstream::parRun(false);
                tmp<GeoField> tsubfld = subsetterPtr().interpolate(fields[i]);
                Pstream::parRun(oldParRun);

                // Send to all processors that don't have a mesh
                for (const int procI : Pstream::subProcs())
                {
                    if (!haveMesh[procI])
                    {
                        OPstream toProc(Pstream::commsTypes::blocking, procI);
                        toProc<< tsubfld();
                    }
                }
            }
        }
    }
    else if (!haveMesh[Pstream::myProcNo()])
    {
        // Don't have mesh (nor fields). Receive empty field from master.

        forAll(masterNames, i)
        {
            const word& name = masterNames[i];

            // Receive field
            IPstream fromMaster
            (
                Pstream::commsTypes::blocking,
                Pstream::masterNo()
            );
            dictionary fieldDict(fromMaster);

            fields.set
            (
                i,
                new GeoField
                (
                    IOobject
                    (
                        name,
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE
                    ),
                    mesh,
                    fieldDict
                )
            );

            //// Write it for next time round (since mesh gets written as well)
            //fields[i].write();
        }
    }
    else
    {
        // Have mesh so just try to load
        forAll(masterNames, i)
        {
            const word& name = masterNames[i];
            IOobject& io = *objects[name];
            io.writeOpt(IOobject::AUTO_WRITE);

            // Load field (but not oldtime)
            readField(io, mesh, i, fields);
        }
    }
}


// Variant of GeometricField::correctBoundaryConditions that only
// evaluates selected patch fields
template<class GeoField, class CoupledPatchType>
void correctCoupledBoundaryConditions(fvMesh& mesh)
{
    HashTable<GeoField*> flds
    (
        mesh.objectRegistry::lookupClass<GeoField>()
    );

    forAllIters(flds, iter)
    {
        GeoField& fld = *iter();

        typename GeoField::Boundary& bfld = fld.boundaryFieldRef();
        if
        (
            Pstream::defaultCommsType == Pstream::commsTypes::blocking
         || Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
        )
        {
            const label nReq = Pstream::nRequests();

            forAll(bfld, patchi)
            {
                auto& pfld = bfld[patchi];
                const auto& fvp = mesh.boundary()[patchi];

                if (fvp.coupled() && !isA<cyclicACMIFvPatch>(fvp))
                {
                    pfld.initEvaluate(Pstream::defaultCommsType);
                }
            }

            // Block for any outstanding requests
            if
            (
                Pstream::parRun()
             && Pstream::defaultCommsType == Pstream::commsTypes::nonBlocking
            )
            {
                Pstream::waitRequests(nReq);
            }

            for (auto& pfld : bfld)
            {
                const auto& fvp = pfld.patch();

                if (fvp.coupled() && !isA<cyclicACMIFvPatch>(fvp))
                {
                    pfld.evaluate(Pstream::defaultCommsType);
                }
            }
        }
        else if (Pstream::defaultCommsType == Pstream::commsTypes::scheduled)
        {
            const lduSchedule& patchSchedule =
                fld.mesh().globalData().patchSchedule();

            forAll(patchSchedule, patchEvali)
            {
                const label patchi = patchSchedule[patchEvali].patch;
                const auto& fvp = mesh.boundary()[patchi];
                auto& pfld = bfld[patchi];

                if (fvp.coupled() && !isA<cyclicACMIFvPatch>(fvp))
                {
                    if (patchSchedule[patchEvali].init)
                    {
                        pfld.initEvaluate(Pstream::commsTypes::scheduled);
                    }
                    else
                    {
                        pfld.evaluate(Pstream::commsTypes::scheduled);
                    }
                }
            }
        }
        else
        {
            FatalErrorInFunction
                << "Unsupported communications type "
                << Pstream::commsTypeNames[Pstream::defaultCommsType]
                << exit(FatalError);
        }
    }
}


// Inplace redistribute mesh and any fields
autoPtr<mapDistributePolyMesh> redistributeAndWrite
(
    autoPtr<fileOperation>&& writeHandler,
    const Time& baseRunTime,
    const boolList& haveMesh,
    const fileName& meshSubDir,
    const bool doReadFields,
    const bool decompose,       // decompose, i.e. read from undecomposed case
    const bool reconstruct,
    const bool overwrite,
    const fileName& proc0CaseName,
    const label nDestProcs,
    const labelList& decomp,
    const fileName& masterInstDir,
    fvMesh& mesh
)
{
    Time& runTime = const_cast<Time&>(mesh.time());
    const bool oldProcCase = runTime.processorCase();

    //// Print some statistics
    //Info<< "Before distribution:" << endl;
    //printMeshData(mesh);


    PtrList<volScalarField> volScalarFields;
    PtrList<volVectorField> volVectorFields;
    PtrList<volSphericalTensorField> volSphereTensorFields;
    PtrList<volSymmTensorField> volSymmTensorFields;
    PtrList<volTensorField> volTensorFields;

    PtrList<surfaceScalarField> surfScalarFields;
    PtrList<surfaceVectorField> surfVectorFields;
    PtrList<surfaceSphericalTensorField> surfSphereTensorFields;
    PtrList<surfaceSymmTensorField> surfSymmTensorFields;
    PtrList<surfaceTensorField> surfTensorFields;

    PtrList<DimensionedField<scalar, volMesh>> dimScalarFields;
    PtrList<DimensionedField<vector, volMesh>> dimVectorFields;
    PtrList<DimensionedField<sphericalTensor, volMesh>> dimSphereTensorFields;
    PtrList<DimensionedField<symmTensor, volMesh>> dimSymmTensorFields;
    PtrList<DimensionedField<tensor, volMesh>> dimTensorFields;

    DynamicList<word> pointFieldNames;


    if (doReadFields)
    {
        // Create 0 sized mesh to do all the generation of zero sized
        // fields on processors that have zero sized meshes. Note that this is
        // only necessary on master but since polyMesh construction with
        // Pstream::parRun does parallel comms we have to do it on all
        // processors
        autoPtr<fvMeshSubset> subsetterPtr;

        const bool allHaveMesh = !haveMesh.found(false);
        if (!allHaveMesh)
        {
            // Find last non-processor patch.
            const polyBoundaryMesh& patches = mesh.boundaryMesh();

            const label nonProcI = (patches.nNonProcessor() - 1);

            if (nonProcI < 0)
            {
                FatalErrorInFunction
                    << "Cannot find non-processor patch on processor "
                    << Pstream::myProcNo() << nl
                    << " Current patches:" << patches.names()
                    << abort(FatalError);
            }

            // Subset 0 cells, no parallel comms.
            // This is used to create zero-sized fields.
            subsetterPtr.reset
            (
                new fvMeshSubset(mesh, bitSet(), nonProcI, false)
            );
        }


        // Get original objects (before incrementing time!)
        if (Pstream::master() && decompose)
        {
            runTime.caseName() = baseRunTime.caseName();
            runTime.processorCase(false);
        }
        IOobjectList objects(mesh, runTime.timeName());
        if (Pstream::master() && decompose)
        {
            runTime.caseName() = proc0CaseName;
            runTime.processorCase(oldProcCase);
        }

        Info<< "From time " << runTime.timeName()
            << " mesh:" << mesh.objectRegistry::objectPath()
            << " have objects:" << objects.names() << endl;

        // We don't want to map the decomposition (mapping already tested when
        // mapping the cell centre field)
        auto iter = objects.find("cellDist");
        if (iter.found())
        {
            objects.erase(iter);
        }


        // volFields

        if (Pstream::master() && decompose)
        {
            runTime.caseName() = baseRunTime.caseName();
            runTime.processorCase(false);
        }
        readFields
        (
            haveMesh,
            mesh,
            subsetterPtr,
            objects,
            volScalarFields
        );

        readFields
        (
            haveMesh,
            mesh,
            subsetterPtr,
            objects,
            volVectorFields
        );

        readFields
        (
            haveMesh,
            mesh,
            subsetterPtr,
            objects,
            volSphereTensorFields
        );

        readFields
        (
            haveMesh,
            mesh,
            subsetterPtr,
            objects,
            volSymmTensorFields
        );

        readFields
        (
            haveMesh,
            mesh,
            subsetterPtr,
            objects,
            volTensorFields
        );


        // surfaceFields

        readFields
        (
            haveMesh,
            mesh,
            subsetterPtr,
            objects,
            surfScalarFields
        );

        readFields
        (
            haveMesh,
            mesh,
            subsetterPtr,
            objects,
            surfVectorFields
        );

        readFields
        (
            haveMesh,
            mesh,
            subsetterPtr,
            objects,
            surfSphereTensorFields
        );

        readFields
        (
            haveMesh,
            mesh,
            subsetterPtr,
            objects,
            surfSymmTensorFields
        );

        readFields
        (
            haveMesh,
            mesh,
            subsetterPtr,
            objects,
            surfTensorFields
        );


        // Dimensioned internal fields
        readFields
        (
            haveMesh,
            mesh,
            subsetterPtr,
            objects,
            dimScalarFields
        );

        readFields
        (
            haveMesh,
            mesh,
            subsetterPtr,
            objects,
            dimVectorFields
        );

        readFields
        (
            haveMesh,
            mesh,
            subsetterPtr,
            objects,
            dimSphereTensorFields
        );

        readFields
        (
            haveMesh,
            mesh,
            subsetterPtr,
            objects,
            dimSymmTensorFields
        );

        readFields
        (
            haveMesh,
            mesh,
            subsetterPtr,
            objects,
            dimTensorFields
        );


        // pointFields currently not supported. Read their names so we
        // can delete them.
        {
            // Get my objects of type
            pointFieldNames.append
            (
                objects.lookupClass(pointScalarField::typeName).sortedNames()
            );
            pointFieldNames.append
            (
                objects.lookupClass(pointVectorField::typeName).sortedNames()
            );
            pointFieldNames.append
            (
                objects.lookupClass
                (
                    pointSphericalTensorField::typeName
                ).sortedNames()
            );
            pointFieldNames.append
            (
                objects.lookupClass
                (
                    pointSymmTensorField::typeName
                ).sortedNames()
            );
            pointFieldNames.append
            (
                objects.lookupClass(pointTensorField::typeName).sortedNames()
            );

            // Make sure all processors have the same set
            Pstream::scatter(pointFieldNames);
        }

        if (Pstream::master() && decompose)
        {
            runTime.caseName() = proc0CaseName;
            runTime.processorCase(oldProcCase);
        }
    }


    // Mesh distribution engine
    fvMeshDistribute distributor(mesh);

    // Do all the distribution of mesh and fields
    autoPtr<mapDistributePolyMesh> rawMap = distributor.distribute(decomp);

    // Print some statistics
    Info<< "After distribution:" << endl;
    printMeshData(mesh);

    // Get other side of processor boundaries
    correctCoupledBoundaryConditions
    <
        volScalarField,
        processorFvPatchField<scalar>
    >(mesh);
    correctCoupledBoundaryConditions
    <
        volVectorField,
        processorFvPatchField<vector>
    >(mesh);
    correctCoupledBoundaryConditions
    <
        volSphericalTensorField,
        processorFvPatchField<sphericalTensor>
    >(mesh);
    correctCoupledBoundaryConditions
    <
        volSymmTensorField,
        processorFvPatchField<symmTensor>
    >(mesh);
    correctCoupledBoundaryConditions
    <
        volTensorField,
        processorFvPatchField<tensor>
    >(mesh);
    // No update surface fields


    // Set the minimum write precision
    IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));


    if (!overwrite)
    {
        ++runTime;
        mesh.setInstance(runTime.timeName());
    }
    else
    {
        mesh.setInstance(masterInstDir);
    }


    IOmapDistributePolyMesh map
    (
        IOobject
        (
            "procAddressing",
            mesh.facesInstance(),
            polyMesh::meshSubDir,
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        )
    );
    map.transfer(rawMap());


    if (reconstruct)
    {
        if (Pstream::master())
        {
            Info<< "Setting caseName to " << baseRunTime.caseName()
                << " to write reconstructed mesh and fields." << endl;
            runTime.caseName() = baseRunTime.caseName();
            const bool oldProcCase(runTime.processorCase(false));
            const bool oldParRun = Pstream::parRun(false);

            mesh.write();
            topoSet::removeFiles(mesh);
            for (const word& fieldName : pointFieldNames)
            {
                IOobject io
                (
                    fieldName,
                    runTime.timeName(),
                    mesh
                );

                const fileName fieldFile(io.objectPath());
                if (topoSet::debug) DebugVar(fieldFile);
                rm(fieldFile);
            }
            Pstream::parRun(oldParRun);

            // Now we've written all. Reset caseName on master
            Info<< "Restoring caseName to " << proc0CaseName << endl;
            runTime.caseName() = proc0CaseName;
            runTime.processorCase(oldProcCase);
        }
    }
    else
    {
        autoPtr<fileOperation> defaultHandler;
        if (writeHandler.valid())
        {
            defaultHandler = fileHandler(std::move(writeHandler));
        }

        mesh.write();

        if (defaultHandler.valid())
        {
            writeHandler = fileHandler(std::move(defaultHandler));
        }
        topoSet::removeFiles(mesh);
        for (const word& fieldName : pointFieldNames)
        {
            IOobject io
            (
                fieldName,
                runTime.timeName(),
                mesh
            );

            const fileName fieldFile(io.objectPath());
            if (topoSet::debug) DebugVar(fieldFile);
            rm(fieldFile);
        }
    }
    Info<< "Written redistributed mesh to " << mesh.facesInstance() << nl
        << endl;


    if (decompose || reconstruct)
    {
        // Decompose (1 -> N) or reconstruct (N -> 1)
        // so {boundary,cell,face,point}ProcAddressing have meaning
        writeProcAddressing(std::move(writeHandler), mesh, map, decompose);
    }
    else
    {
        // Redistribute (N -> M)
        // {boundary,cell,face,point}ProcAddressing would be incorrect
        // - can either remove or redistribute previous
        removeProcAddressing(mesh);
    }


    // Refinement data
    {

        // Read refinement data
        if (Pstream::master() && decompose)
        {
            runTime.caseName() = baseRunTime.caseName();
            runTime.processorCase(false);
        }
        IOobject io
        (
            "dummy",
            mesh.facesInstance(),
            polyMesh::meshSubDir,
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            false
        );

        hexRef8Data refData(io);
        if (Pstream::master() && decompose)
        {
            runTime.caseName() = proc0CaseName;
            runTime.processorCase(oldProcCase);
        }

        // Make sure all processors have valid data (since only some will
        // read)
        refData.sync(io);

        // Distribute
        refData.distribute(map);


        // Now we've read refinement data we can remove it
        meshRefinement::removeFiles(mesh);

        if (reconstruct)
        {
            if (Pstream::master())
            {
                const bool oldParRun = Pstream::parRun(false);

                Info<< "Setting caseName to " << baseRunTime.caseName()
                    << " to write reconstructed refinement data." << endl;
                runTime.caseName() = baseRunTime.caseName();
                const bool oldProcCase(runTime.processorCase(false));

                refData.write();

                // Now we've written all. Reset caseName on master
                Info<< "Restoring caseName to " << proc0CaseName << endl;
                runTime.caseName() = proc0CaseName;
                runTime.processorCase(oldProcCase);

                Pstream::parRun(oldParRun);
            }
        }
        else
        {
            refData.write();
        }
    }

    //// Sets. Disabled for now.
    //{
    //    // Read sets
    //    if (Pstream::master() && decompose)
    //    {
    //        runTime.caseName() = baseRunTime.caseName();
    //        runTime.processorCase(false);
    //    }
    //    IOobjectList objects(mesh, mesh.facesInstance(), "polyMesh/sets");
    //
    //    PtrList<cellSet> cellSets;
    //    ReadFields(objects, cellSets);
    //
    //    if (Pstream::master() && decompose)
    //    {
    //        runTime.caseName() = proc0CaseName;
    //        runTime.processorCase(oldProcCase);
    //    }
    //
    //    forAll(cellSets, i)
    //    {
    //        cellSets[i].distribute(map);
    //    }
    //
    //    if (reconstruct)
    //    {
    //        if (Pstream::master())
    //        {
    //            Info<< "Setting caseName to " << baseRunTime.caseName()
    //                << " to write reconstructed refinement data." << endl;
    //            runTime.caseName() = baseRunTime.caseName();
    //            const bool oldProcCase(runTime.processorCase(false));
    //
    //            forAll(cellSets, i)
    //            {
    //                cellSets[i].distribute(map);
    //            }
    //
    //            // Now we've written all. Reset caseName on master
    //            Info<< "Restoring caseName to " << proc0CaseName << endl;
    //            runTime.caseName() = proc0CaseName;
    //            runTime.processorCase(oldProcCase);
    //        }
    //    }
    //    else
    //    {
    //        forAll(cellSets, i)
    //        {
    //            cellSets[i].distribute(map);
    //        }
    //    }
    //}


    return autoPtr<mapDistributePolyMesh>::New(std::move(map));
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
// Field Mapping
//
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

autoPtr<mapDistributePolyMesh> createReconstructMap
(
    const autoPtr<fvMesh>& baseMeshPtr,
    const fvMesh& mesh,
    const labelList& cellProcAddressing,
    const labelList& faceProcAddressing,
    const labelList& pointProcAddressing,
    const labelList& boundaryProcAddressing
)
{
    // Send addressing to master
    labelListList cellAddressing(Pstream::nProcs());
    cellAddressing[Pstream::myProcNo()] = cellProcAddressing;
    Pstream::gatherList(cellAddressing);

    labelListList faceAddressing(Pstream::nProcs());
    faceAddressing[Pstream::myProcNo()] = faceProcAddressing;
    Pstream::gatherList(faceAddressing);

    labelListList pointAddressing(Pstream::nProcs());
    pointAddressing[Pstream::myProcNo()] = pointProcAddressing;
    Pstream::gatherList(pointAddressing);

    labelListList boundaryAddressing(Pstream::nProcs());
    {
        // Remove -1 entries
        DynamicList<label> patchProcAddressing(boundaryProcAddressing.size());
        forAll(boundaryProcAddressing, i)
        {
            if (boundaryProcAddressing[i] != -1)
            {
                patchProcAddressing.append(boundaryProcAddressing[i]);
            }
        }
        boundaryAddressing[Pstream::myProcNo()] = patchProcAddressing;
        Pstream::gatherList(boundaryAddressing);
    }


    autoPtr<mapDistributePolyMesh> mapPtr;

    if (baseMeshPtr && baseMeshPtr->nCells())
    {
        const fvMesh& baseMesh = *baseMeshPtr;

        labelListList cellSubMap(Pstream::nProcs());
        cellSubMap[Pstream::masterNo()] = identity(mesh.nCells());

        mapDistribute cellMap
        (
            baseMesh.nCells(),
            std::move(cellSubMap),
            std::move(cellAddressing)
        );

        labelListList faceSubMap(Pstream::nProcs());
        faceSubMap[Pstream::masterNo()] = identity(mesh.nFaces());

        mapDistribute faceMap
        (
            baseMesh.nFaces(),
            std::move(faceSubMap),
            std::move(faceAddressing),
            false,          //subHasFlip
            true            //constructHasFlip
        );

        labelListList pointSubMap(Pstream::nProcs());
        pointSubMap[Pstream::masterNo()] = identity(mesh.nPoints());

        mapDistribute pointMap
        (
            baseMesh.nPoints(),
            std::move(pointSubMap),
            std::move(pointAddressing)
        );

        labelListList patchSubMap(Pstream::nProcs());
        // Send (filtered) patches to master
        patchSubMap[Pstream::masterNo()] =
            boundaryAddressing[Pstream::myProcNo()];

        mapDistribute patchMap
        (
            baseMesh.boundaryMesh().size(),
            std::move(patchSubMap),
            std::move(boundaryAddressing)
        );

        const label nOldPoints = mesh.nPoints();
        const label nOldFaces = mesh.nFaces();
        const label nOldCells = mesh.nCells();

        const polyBoundaryMesh& pbm = mesh.boundaryMesh();
        labelList oldPatchStarts(pbm.size());
        labelList oldPatchNMeshPoints(pbm.size());
        forAll(pbm, patchI)
        {
            oldPatchStarts[patchI] = pbm[patchI].start();
            oldPatchNMeshPoints[patchI] = pbm[patchI].nPoints();
        }

        mapPtr.reset
        (
            new mapDistributePolyMesh
            (
                nOldPoints,
                nOldFaces,
                nOldCells,
                std::move(oldPatchStarts),
                std::move(oldPatchNMeshPoints),
                std::move(pointMap),
                std::move(faceMap),
                std::move(cellMap),
                std::move(patchMap)
            )
        );
    }
    else
    {
        labelListList cellSubMap(Pstream::nProcs());
        cellSubMap[Pstream::masterNo()] = identity(mesh.nCells());
        labelListList cellConstructMap(Pstream::nProcs());

        mapDistribute cellMap
        (
            0,
            std::move(cellSubMap),
            std::move(cellConstructMap)
        );

        labelListList faceSubMap(Pstream::nProcs());
        faceSubMap[Pstream::masterNo()] = identity(mesh.nFaces());
        labelListList faceConstructMap(Pstream::nProcs());

        mapDistribute faceMap
        (
            0,
            std::move(faceSubMap),
            std::move(faceConstructMap),
            false,          //subHasFlip
            true            //constructHasFlip
        );

        labelListList pointSubMap(Pstream::nProcs());
        pointSubMap[Pstream::masterNo()] = identity(mesh.nPoints());
        labelListList pointConstructMap(Pstream::nProcs());

        mapDistribute pointMap
        (
            0,
            std::move(pointSubMap),
            std::move(pointConstructMap)
        );

        labelListList patchSubMap(Pstream::nProcs());
        // Send (filtered) patches to master
        patchSubMap[Pstream::masterNo()] =
            boundaryAddressing[Pstream::myProcNo()];
        labelListList patchConstructMap(Pstream::nProcs());

        mapDistribute patchMap
        (
            0,
            std::move(patchSubMap),
            std::move(patchConstructMap)
        );

        const label nOldPoints = mesh.nPoints();
        const label nOldFaces = mesh.nFaces();
        const label nOldCells = mesh.nCells();

        const polyBoundaryMesh& pbm = mesh.boundaryMesh();
        labelList oldPatchStarts(pbm.size());
        labelList oldPatchNMeshPoints(pbm.size());
        forAll(pbm, patchI)
        {
            oldPatchStarts[patchI] = pbm[patchI].start();
            oldPatchNMeshPoints[patchI] = pbm[patchI].nPoints();
        }

        mapPtr.reset
        (
            new mapDistributePolyMesh
            (
                nOldPoints,
                nOldFaces,
                nOldCells,
                std::move(oldPatchStarts),
                std::move(oldPatchNMeshPoints),
                std::move(pointMap),
                std::move(faceMap),
                std::move(cellMap),
                std::move(patchMap)
            )
        );
    }

    return mapPtr;
}


void readProcAddressing
(
    const fvMesh& mesh,
    const autoPtr<fvMesh>& baseMeshPtr,
    autoPtr<mapDistributePolyMesh>& distMap
)
{
    //IOobject io
    //(
    //    "procAddressing",
    //    mesh.facesInstance(),
    //    polyMesh::meshSubDir,
    //    mesh,
    //    IOobject::MUST_READ
    //);
    //if (io.typeHeaderOk<labelIOList>(true))
    //{
    //    Pout<< "Reading addressing from " << io.name() << " at "
    //        << mesh.facesInstance() << nl << endl;
    //    distMap.reset(new IOmapDistributePolyMesh(io));
    //}
    //else
    {
        Info<< "Reading addressing from procXXXAddressing at "
            << mesh.facesInstance() << nl << endl;
        labelIOList cellProcAddressing
        (
            IOobject
            (
                "cellProcAddressing",
                mesh.facesInstance(),
                polyMesh::meshSubDir,
                mesh,
                IOobject::READ_IF_PRESENT
            ),
            labelList()
        );
        labelIOList faceProcAddressing
        (
            IOobject
            (
                "faceProcAddressing",
                mesh.facesInstance(),
                polyMesh::meshSubDir,
                mesh,
                IOobject::READ_IF_PRESENT
            ),
            labelList()
        );
        labelIOList pointProcAddressing
        (
            IOobject
            (
                "pointProcAddressing",
                mesh.facesInstance(),
                polyMesh::meshSubDir,
                mesh,
                IOobject::READ_IF_PRESENT
            ),
            labelList()
        );
        labelIOList boundaryProcAddressing
        (
            IOobject
            (
                "boundaryProcAddressing",
                mesh.facesInstance(),
                polyMesh::meshSubDir,
                mesh,
                IOobject::READ_IF_PRESENT
            ),
            labelList()
        );


        if
        (
            mesh.nCells() != cellProcAddressing.size()
         || mesh.nPoints() != pointProcAddressing.size()
         || mesh.nFaces() != faceProcAddressing.size()
         || mesh.boundaryMesh().size() != boundaryProcAddressing.size()
        )
        {
            FatalErrorInFunction
                << "Read addressing inconsistent with mesh sizes" << nl
                << "cells:" << mesh.nCells()
                << " addressing:" << cellProcAddressing.objectPath()
                << " size:" << cellProcAddressing.size() << nl
                << "faces:" << mesh.nFaces()
                << " addressing:" << faceProcAddressing.objectPath()
                << " size:" << faceProcAddressing.size() << nl
                << "points:" << mesh.nPoints()
                << " addressing:" << pointProcAddressing.objectPath()
                << " size:" << pointProcAddressing.size()
                << "patches:" << mesh.boundaryMesh().size()
                << " addressing:" << boundaryProcAddressing.objectPath()
                << " size:" << boundaryProcAddressing.size()
                << exit(FatalError);
        }

        distMap.clear();
        distMap = createReconstructMap
        (
            baseMeshPtr,
            mesh,
            cellProcAddressing,
            faceProcAddressing,
            pointProcAddressing,
            boundaryProcAddressing
        );
    }
}


void reconstructMeshFields
(
    const parFvFieldReconstructor& fvReconstructor,
    const IOobjectList& objects,
    const wordRes& selectedFields
)
{
    // Dimensioned fields

    fvReconstructor.reconstructFvVolumeInternalFields<scalar>
    (
        objects,
        selectedFields
    );
    fvReconstructor.reconstructFvVolumeInternalFields<vector>
    (
        objects,
        selectedFields
    );
    fvReconstructor.reconstructFvVolumeInternalFields<sphericalTensor>
    (
        objects,
        selectedFields
    );
    fvReconstructor.reconstructFvVolumeInternalFields<symmTensor>
    (
        objects,
        selectedFields
    );
    fvReconstructor.reconstructFvVolumeInternalFields<tensor>
    (
        objects,
        selectedFields
    );


    // volFields

    fvReconstructor.reconstructFvVolumeFields<scalar>
    (
        objects,
        selectedFields
    );
    fvReconstructor.reconstructFvVolumeFields<vector>
    (
        objects,
        selectedFields
    );
    fvReconstructor.reconstructFvVolumeFields<sphericalTensor>
    (
        objects,
        selectedFields
    );
    fvReconstructor.reconstructFvVolumeFields<symmTensor>
    (
        objects,
        selectedFields
    );
    fvReconstructor.reconstructFvVolumeFields<tensor>
    (
        objects,
        selectedFields
    );


    // surfaceFields

    fvReconstructor.reconstructFvSurfaceFields<scalar>
    (
        objects,
        selectedFields
    );
    fvReconstructor.reconstructFvSurfaceFields<vector>
    (
        objects,
        selectedFields
    );
    fvReconstructor.reconstructFvSurfaceFields<sphericalTensor>
    (
        objects,
        selectedFields
    );
    fvReconstructor.reconstructFvSurfaceFields<symmTensor>
    (
        objects,
        selectedFields
    );
    fvReconstructor.reconstructFvSurfaceFields<tensor>
    (
        objects,
        selectedFields
    );
}


void reconstructLagrangian
(
    autoPtr<parLagrangianRedistributor>& lagrangianReconstructorPtr,
    const fvMesh& baseMesh,
    const fvMesh& mesh,
    const mapDistributePolyMesh& distMap,
    const wordRes& selectedLagrangianFields
)
{
    // Clouds (note: might not be present on all processors)

    wordList cloudNames;
    List<wordList> fieldNames;
    // Find all cloudNames on all processors
    parLagrangianRedistributor::findClouds(mesh, cloudNames, fieldNames);

    if (cloudNames.size())
    {
        if (!lagrangianReconstructorPtr)
        {
            lagrangianReconstructorPtr.reset
            (
                new parLagrangianRedistributor
                (
                    mesh,
                    baseMesh,
                    mesh.nCells(),      // range of cell indices in clouds
                    distMap
                )
            );
        }
        const parLagrangianRedistributor& lagrangianReconstructor =
            *lagrangianReconstructorPtr;

        for (const word& cloudName : cloudNames)
        {
            Info<< "Reconstructing lagrangian fields for cloud "
                << cloudName << nl << endl;

            autoPtr<mapDistributeBase> lagrangianMapPtr =
                lagrangianReconstructor.redistributeLagrangianPositions
                (
                    cloudName
                );

            const mapDistributeBase& lagrangianMap = *lagrangianMapPtr;

            IOobjectList cloudObjs
            (
                mesh,
                mesh.time().timeName(),
                cloud::prefix/cloudName
            );

            lagrangianReconstructor.redistributeFields<label>
            (
                lagrangianMap,
                cloudName,
                cloudObjs,
                selectedLagrangianFields
            );
            lagrangianReconstructor.redistributeFieldFields<label>
            (
                lagrangianMap,
                cloudName,
                cloudObjs,
                selectedLagrangianFields
            );
            lagrangianReconstructor.redistributeFields<scalar>
            (
                lagrangianMap,
                cloudName,
                cloudObjs,
                selectedLagrangianFields
            );
            lagrangianReconstructor.redistributeFieldFields<scalar>
            (
                lagrangianMap,
                cloudName,
                cloudObjs,
                selectedLagrangianFields
            );
            lagrangianReconstructor.redistributeFields<vector>
            (
                lagrangianMap,
                cloudName,
                cloudObjs,
                selectedLagrangianFields
            );
            lagrangianReconstructor.redistributeFieldFields<vector>
            (
                lagrangianMap,
                cloudName,
                cloudObjs,
                selectedLagrangianFields
            );
            lagrangianReconstructor.redistributeFields
            <sphericalTensor>
            (
                lagrangianMap,
                cloudName,
                cloudObjs,
                selectedLagrangianFields
            );
            lagrangianReconstructor.redistributeFieldFields
            <sphericalTensor>
            (
                lagrangianMap,
                cloudName,
                cloudObjs,
                selectedLagrangianFields
            );
            lagrangianReconstructor.redistributeFields<symmTensor>
            (
                lagrangianMap,
                cloudName,
                cloudObjs,
                selectedLagrangianFields
            );
            lagrangianReconstructor.redistributeFieldFields
            <symmTensor>
            (
                lagrangianMap,
                cloudName,
                cloudObjs,
                selectedLagrangianFields
            );
            lagrangianReconstructor.redistributeFields<tensor>
            (
                lagrangianMap,
                cloudName,
                cloudObjs,
                selectedLagrangianFields
            );
            lagrangianReconstructor.redistributeFieldFields<tensor>
            (
                lagrangianMap,
                cloudName,
                cloudObjs,
                selectedLagrangianFields
            );
        }
    }
}


// Read clouds (note: might not be present on all processors)
void readLagrangian
(
    const fvMesh& mesh,
    const wordList& cloudNames,
    const wordRes& selectedLagrangianFields,
    PtrList<unmappedPassivePositionParticleCloud>& clouds
)
{
    (void)mesh.tetBasePtIs();

    forAll(cloudNames, i)
    {
        //Pout<< "Loading cloud " << cloudNames[i] << endl;
        clouds.set
        (
            i,
            new unmappedPassivePositionParticleCloud(mesh, cloudNames[i], false)
        );


        //for (passivePositionParticle& p : clouds[i]))
        //{
        //    Pout<< "Particle position:" << p.position()
        //        << " cell:" << p.cell()
        //        << " with cc:" << mesh.cellCentres()[p.cell()]
        //        << endl;
        //}


        IOobjectList cloudObjs(clouds[i], clouds[i].time().timeName());

        //Pout<< "Found clould objects:" << cloudObjs.names() << endl;

        parLagrangianRedistributor::readFields
        <IOField<label>>
        (
            clouds[i],
            cloudObjs,
            selectedLagrangianFields
        );
        parLagrangianRedistributor::readFields
        <IOField<Field<label>>>
        (
            clouds[i],
            cloudObjs,
            selectedLagrangianFields
        );
        parLagrangianRedistributor::readFields
        <CompactIOField<Field<label>, label>>
        (
            clouds[i],
            cloudObjs,
            selectedLagrangianFields
        );


        parLagrangianRedistributor::readFields
        <IOField<scalar>>
        (
            clouds[i],
            cloudObjs,
            selectedLagrangianFields
        );
        parLagrangianRedistributor::readFields
        <IOField<Field<scalar>>>
        (
            clouds[i],
            cloudObjs,
            selectedLagrangianFields
        );
        parLagrangianRedistributor::readFields
        <CompactIOField<Field<scalar>, scalar>>
        (
            clouds[i],
            cloudObjs,
            selectedLagrangianFields
        );


        parLagrangianRedistributor::readFields
        <IOField<vector>>
        (
            clouds[i],
            cloudObjs,
            selectedLagrangianFields
        );
        parLagrangianRedistributor::readFields
        <IOField<Field<vector>>>
        (
            clouds[i],
            cloudObjs,
            selectedLagrangianFields
        );
        parLagrangianRedistributor::readFields
        <CompactIOField<Field<vector>, vector>>
        (
            clouds[i],
            cloudObjs,
            selectedLagrangianFields
        );


        parLagrangianRedistributor::readFields
        <IOField<sphericalTensor>>
        (
            clouds[i],
            cloudObjs,
            selectedLagrangianFields
        );
        parLagrangianRedistributor::readFields
        <IOField<Field<sphericalTensor>>>
        (
            clouds[i],
            cloudObjs,
            selectedLagrangianFields
        );
        parLagrangianRedistributor::readFields
        <CompactIOField<Field<sphericalTensor>, sphericalTensor>>
        (
            clouds[i],
            cloudObjs,
            selectedLagrangianFields
        );


        parLagrangianRedistributor::readFields
        <IOField<symmTensor>>
        (
            clouds[i],
            cloudObjs,
            selectedLagrangianFields
        );
        parLagrangianRedistributor::readFields
        <IOField<Field<symmTensor>>>
        (
            clouds[i],
            cloudObjs,
            selectedLagrangianFields
        );
        parLagrangianRedistributor::readFields
        <CompactIOField<Field<symmTensor>, symmTensor>>
        (
            clouds[i],
            cloudObjs,
            selectedLagrangianFields
        );


        parLagrangianRedistributor::readFields
        <IOField<tensor>>
        (
            clouds[i],
            cloudObjs,
            selectedLagrangianFields
        );
        parLagrangianRedistributor::readFields
        <IOField<Field<tensor>>>
        (
            clouds[i],
            cloudObjs,
            selectedLagrangianFields
        );
        parLagrangianRedistributor::readFields
        <CompactIOField<Field<tensor>, tensor>>
        (
            clouds[i],
            cloudObjs,
            selectedLagrangianFields
        );
    }
}


void redistributeLagrangian
(
    autoPtr<parLagrangianRedistributor>& lagrangianReconstructorPtr,
    const fvMesh& mesh,
    const label nOldCells,
    const mapDistributePolyMesh& distMap,
    PtrList<unmappedPassivePositionParticleCloud>& clouds
)
{
    if (clouds.size())
    {
        if (!lagrangianReconstructorPtr)
        {
            lagrangianReconstructorPtr.reset
            (
                new parLagrangianRedistributor
                (
                    mesh,
                    mesh,
                    nOldCells,  // range of cell indices in clouds
                    distMap
                )
            );
        }
        const parLagrangianRedistributor& distributor =
            lagrangianReconstructorPtr();

        forAll(clouds, i)
        {
            autoPtr<mapDistributeBase> lagrangianMapPtr =
                distributor.redistributeLagrangianPositions(clouds[i]);
            const mapDistributeBase& lagrangianMap = *lagrangianMapPtr;

            distributor.redistributeStoredFields
            <IOField<label>>
            (
                lagrangianMap,
                clouds[i]
            );
            distributor.redistributeStoredFields
            <IOField<Field<label>>>
            (
                lagrangianMap,
                clouds[i]
            );
            distributor.redistributeStoredFields
            <CompactIOField<Field<label>, label>>
            (
                lagrangianMap,
                clouds[i]
            );


            distributor.redistributeStoredFields
            <IOField<scalar>>
            (
                lagrangianMap,
                clouds[i]
            );
            distributor.redistributeStoredFields
            <IOField<Field<scalar>>>
            (
                lagrangianMap,
                clouds[i]
            );
            distributor.redistributeStoredFields
            <CompactIOField<Field<scalar>, scalar>>
            (
                lagrangianMap,
                clouds[i]
            );


            distributor.redistributeStoredFields
            <IOField<vector>>
            (
                lagrangianMap,
                clouds[i]
            );
            distributor.redistributeStoredFields
            <IOField<Field<vector>>>
            (
                lagrangianMap,
                clouds[i]
            );
            distributor.redistributeStoredFields
            <CompactIOField<Field<vector>, vector>>
            (
                lagrangianMap,
                clouds[i]
            );


            distributor.redistributeStoredFields
            <IOField<sphericalTensor>>
            (
                lagrangianMap,
                clouds[i]
            );
            distributor.redistributeStoredFields
            <IOField<Field<sphericalTensor>>>
            (
                lagrangianMap,
                clouds[i]
            );
            distributor.redistributeStoredFields
            <CompactIOField<Field<sphericalTensor>, sphericalTensor>>
            (
                lagrangianMap,
                clouds[i]
            );


            distributor.redistributeStoredFields
            <IOField<symmTensor>>
            (
                lagrangianMap,
                clouds[i]
            );
            distributor.redistributeStoredFields
            <IOField<Field<symmTensor>>>
            (
                lagrangianMap,
                clouds[i]
            );
            distributor.redistributeStoredFields
            <CompactIOField<Field<symmTensor>, symmTensor>>
            (
                lagrangianMap,
                clouds[i]
            );


            distributor.redistributeStoredFields
            <IOField<tensor>>
            (
                lagrangianMap,
                clouds[i]
            );
            distributor.redistributeStoredFields
            <IOField<Field<tensor>>>
            (
                lagrangianMap,
                clouds[i]
            );
            distributor.redistributeStoredFields
            <CompactIOField<Field<tensor>, tensor>>
            (
                lagrangianMap,
                clouds[i]
            );
        }
    }
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Redistribute decomposed mesh and fields according"
        " to the decomposeParDict settings.\n"
        "Optionally run in decompose/reconstruct mode"
    );

    argList::noFunctionObjects();  // Never use function objects

    // enable -constant ... if someone really wants it
    // enable -zeroTime to prevent accidentally trashing the initial fields
    timeSelector::addOptions(true, true);

    #include "addAllRegionOptions.H"

    #include "addOverwriteOption.H"
    argList::addBoolOption("decompose", "Decompose case");
    argList::addBoolOption("reconstruct", "Reconstruct case");
    argList::addDryRunOption
    (
        "Test without writing the decomposition. "
        "Changes -cellDist to only write volScalarField."
    );
    argList::addBoolOption
    (
        "cellDist",
        "Write cell distribution as a labelList - for use with 'manual' "
        "decomposition method or as a volScalarField for post-processing."
    );
    argList::addBoolOption
    (
        "newTimes",
        "Only reconstruct new times (i.e. that do not exist already)"
    );


    // Handle arguments
    // ~~~~~~~~~~~~~~~~
    // (replacement for setRootCase that does not abort)

    argList args(argc, argv);


    // As much as possible avoid synchronised operation. To be looked at more
    // closely for the three scenarios:
    // - decompose - reads on master (and from parent directory) and sends
    //               dictionary to slaves
    // - distribute - reads on potentially a different number of processors
    //                than it writes to
    // - reconstruct - reads parallel, write on master only and to parent
    //                 directory
    autoPtr<fileOperation> writeHandler;
    if
    (
        fileHandler().type()
     != fileOperations::uncollatedFileOperation::typeName
    )
    {
        // Install 'uncollated' as fileHandler. Save old one in writeHandler.
        writeHandler = fileHandler(fileOperation::NewUncollated());
    }

    // Switch off parallel synchronisation of cached time directories
    fileHandler().distributed(true);

    // File handler to be used for writing
    const fileOperation& fh
    (
        writeHandler.valid()
      ? writeHandler()
      : fileHandler()
    );



    // Make sure to call findTimes on all processors to force caching of
    // time directories
    (void)fileHandler().findTimes(args.path(), "constant");


    #include "foamDlOpenLibs.H"

    const bool reconstruct = args.found("reconstruct");
    const bool writeCellDist = args.found("cellDist");
    const bool dryrun = args.dryRun();
    const bool newTimes = args.found("newTimes");

    bool decompose = args.found("decompose");
    bool overwrite = args.found("overwrite");

    // Disable NaN setting and floating point error trapping. This is to avoid
    // any issues inside the field redistribution inside fvMeshDistribute
    // which temporarily moves processor faces into existing patches. These
    // will now not have correct values. After all bits have been assembled
    // the processor fields will be restored but until then there might
    // be uninitialised values which might upset some patch field constructors.
    // Workaround by disabling floating point error trapping. TBD: have
    // actual field redistribution instead of subsetting inside
    // fvMeshDistribute.
    Foam::sigFpe::unset(true);

    const wordRes selectedFields;
    const wordRes selectedLagrangianFields;


    if (decompose)
    {
        Info<< "Decomposing case (like decomposePar)" << nl << endl;
        if (reconstruct)
        {
            FatalErrorInFunction
                << "Cannot specify both -decompose and -reconstruct"
                << exit(FatalError);
        }
    }
    else if (reconstruct)
    {
        Info<< "Reconstructing case (like reconstructParMesh)" << nl << endl;
    }


    if (decompose || reconstruct)
    {
        if (!overwrite)
        {
            WarningInFunction
                << "Working in decompose or reconstruction mode automatically"
                << " implies -overwrite" << nl << endl;
            overwrite = true;
        }
    }


    if (!Pstream::parRun())
    {
        FatalErrorInFunction
            << ": This utility can only be run parallel"
            << exit(FatalError);
    }


    if (!isDir(args.rootPath()))
    {
        FatalErrorInFunction
            << ": cannot open root directory " << args.rootPath()
            << exit(FatalError);
    }

    // Detect if running data-distributed (multiple roots)
    bool nfs = true;
    {
        List<fileName> roots(1, args.rootPath());
        combineReduce(roots, ListOps::uniqueEqOp<fileName>());
        nfs = (roots.size() == 1);
    }

    if (!nfs)
    {
        Info<< "Detected multiple roots i.e. non-nfs running"
            << nl << endl;
    }

    // Check if we have processor directories. Ideally would like to
    // use fileHandler().dirPath here but we don't have runTime yet and
    // want to delay constructing runTime until we've synced all time
    // directories...
    const fileName procDir(fileHandler().filePath(args.path()));
    if (isDir(procDir))
    {
        if (decompose)
        {
            Info<< "Removing existing processor directory" << procDir << endl;
            fileHandler().rmDir(procDir);
        }
    }
    else
    {
        // Directory does not exist. If this happens on master -> decompose mode
        if (Pstream::master())
        {
           decompose = true;
           Info<< "No processor directories; switching on decompose mode"
               << nl << endl;
        }
    }
    // If master changed to decompose mode make sure all nodes know about it
    Pstream::scatter(decompose);


    // If running distributed we have problem of new processors not finding
    // a system/controlDict. However if we switch on the master-only reading
    // the problem becomes that the time directories are differing sizes and
    // e.g. latestTime will pick up a different time (which causes createTime.H
    // to abort). So for now make sure to have master times on all
    // processors
    if (!reconstruct)
    {
        Info<< "Creating time directories on all processors" << nl << endl;
        createTimeDirs(args.path());
    }

    // Construct time
    // ~~~~~~~~~~~~~~

    #include "createTime.H"
    runTime.functionObjects().off();  // Extra safety?


    // Save local processor0 casename
    const fileName proc0CaseName = runTime.caseName();
    const bool oldProcCase = runTime.processorCase();


    // Construct undecomposed Time
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // This will read the same controlDict but might have a different
    // set of times so enforce same times

    if (!nfs)
    {
        Info<< "Creating time directories for undecomposed Time"
            << " on all processors" << nl << endl;
        createTimeDirs(args.globalPath());
    }


    Info<< "Create undecomposed database"<< nl << endl;
    Time baseRunTime
    (
        runTime.controlDict(),
        runTime.rootPath(),
        runTime.globalCaseName(),
        runTime.system(),
        runTime.constant(),
        false                   // enableFunctionObjects
    );


    wordHashSet masterTimeDirSet;
    if (newTimes)
    {
        instantList baseTimeDirs(baseRunTime.times());
        for (const instant& t : baseTimeDirs)
        {
            masterTimeDirSet.insert(t.name());
        }
    }


    // Allow override of decomposeParDict location
    const fileName decompDictFile =
        args.getOrDefault<fileName>("decomposeParDict", "");

    // Get region names
    #include "getAllRegionOptions.H"

    if (regionNames.size() == 1 && regionNames[0] != polyMesh::defaultRegion)
    {
        Info<< "Using region: " << regionNames[0] << nl << endl;
    }


    // Demand driven lagrangian mapper
    autoPtr<parLagrangianRedistributor> lagrangianReconstructorPtr;


    if (reconstruct)
    {
        // use the times list from the master processor
        // and select a subset based on the command-line options
        instantList timeDirs = timeSelector::select(runTime.times(), args);
        Pstream::scatter(timeDirs);

        if (timeDirs.empty())
        {
            FatalErrorInFunction
                << "No times selected"
                << exit(FatalError);
        }


        // Pass1 : reconstruct mesh and addressing
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


        Info<< nl
            << "Pass1 : reconstructing mesh and addressing" << nl << endl;


        forAll(regionNames, regioni)
        {
            const word& regionName = regionNames[regioni];
            const word& regionDir =
            (
                regionName == polyMesh::defaultRegion ? word::null : regionName
            );
            const fileName meshSubDir(regionDir/polyMesh::meshSubDir);

            Info<< "\n\nReconstructing mesh " << regionName << nl << endl;

            // Loop over all times
            forAll(timeDirs, timeI)
            {
                // Set time for global database
                runTime.setTime(timeDirs[timeI], timeI);
                baseRunTime.setTime(timeDirs[timeI], timeI);

                Info<< "Time = " << runTime.timeName() << endl << endl;


                // See where the mesh is
                fileName facesInstance = runTime.findInstance
                (
                    meshSubDir,
                    "faces",
                    IOobject::READ_IF_PRESENT
                );
                //Pout<< "facesInstance:" << facesInstance << endl;

                Pstream::scatter(facesInstance);

                // Check who has a mesh (by checking for 'faces' file)
                const boolList haveMesh
                (
                    haveFacesFile
                    (
                        runTime.path()/facesInstance/meshSubDir
                    )
                );


                // Addressing back to reconstructed mesh as xxxProcAddressing.
                // - all processors have consistent faceProcAddressing
                // - processors with no mesh don't need faceProcAddressing


                // Note: filePath searches up on processors that don't have
                //       processor if instance = constant so explicitly check
                //       found filename.
                bool haveAddressing = false;
                if (haveMesh[Pstream::myProcNo()])
                {
                    // Read faces (just to know their size)
                    faceCompactIOList faces
                    (
                        IOobject
                        (
                            "faces",
                            facesInstance,
                            meshSubDir,
                            runTime,
                            IOobject::MUST_READ
                        )
                    );

                    // Check faceProcAddressing
                    labelIOList faceProcAddressing
                    (
                        IOobject
                        (
                            "faceProcAddressing",
                            facesInstance,
                            meshSubDir,
                            runTime,
                            IOobject::READ_IF_PRESENT
                        ),
                        labelList()
                    );
                    if
                    (
                        faceProcAddressing.headerOk()
                     && faceProcAddressing.size() == faces.size()
                    )
                    {
                        haveAddressing = true;
                    }
                }
                else
                {
                    // Have no mesh. Don't need addressing
                    haveAddressing = true;
                }


                // Additionally check for master faces being readable. Could
                // do even more checks, e.g. global number of cells same
                // as cellProcAddressing
                bool haveUndecomposedMesh = false;
                if (Pstream::master())
                {
                    Info<< "Checking " << baseRunTime.caseName()
                        << " for undecomposed mesh" << endl;

                    const bool oldParRun = Pstream::parRun(false);
                    faceCompactIOList facesIO
                    (
                        IOobject
                        (
                            "faces",
                            facesInstance,
                            meshSubDir,
                            baseRunTime,
                            IOobject::NO_READ
                        )
                    );
                    haveUndecomposedMesh = facesIO.headerOk();
                    Pstream::parRun(oldParRun);
                }
                Pstream::scatter(haveUndecomposedMesh);


                if
                (
                    !haveUndecomposedMesh
                 || !returnReduce(haveAddressing, andOp<bool>())
                )
                {
                    Info<< "loading mesh from " << facesInstance << endl;
                    autoPtr<fvMesh> meshPtr = loadOrCreateMesh
                    (
                        decompose,
                        IOobject
                        (
                            regionName,
                            facesInstance,
                            runTime,
                            Foam::IOobject::MUST_READ
                        )
                    );
                    fvMesh& mesh = meshPtr();

                    // Use basic geometry calculation to avoid synchronisation
                    // problems. See comment in routine
                    setBasicGeometry(mesh);


                    // Determine decomposition
                    // ~~~~~~~~~~~~~~~~~~~~~~~

                    Info<< "Reconstructing mesh for time " << facesInstance
                        << endl;

                    label nDestProcs = 1;
                    labelList finalDecomp = labelList(mesh.nCells(), Zero);

                    redistributeAndWrite
                    (
                        std::move(writeHandler),
                        baseRunTime,
                        haveMesh,
                        meshSubDir,
                        false,      // do not read fields
                        false,      // do not read undecomposed case on proc0
                        true,       // write redistributed files to proc0
                        overwrite,
                        proc0CaseName,
                        nDestProcs,
                        finalDecomp,
                        facesInstance,
                        mesh
                    );
                }
            }

            // Make sure all is finished writing until re-reading in pass2
            // below
            fileHandler().flush();


            // Pass2 : read mesh and addressing and reconstruct fields
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            Info<< nl
                << "Pass2 : reconstructing fields" << nl << endl;

            runTime.setTime(timeDirs[0], 0);
            baseRunTime.setTime(timeDirs[0], 0);
            Info<< "Time = " << runTime.timeName() << endl << endl;


            // Read undecomposed mesh on master and 'empty' mesh
            // (zero faces, point, cells but valid patches and zones) on slaves.
            // This is a bit of tricky code and hidden inside fvMeshTools for
            // now.
            Info<< "Reading undecomposed mesh (on master)" << endl;
            autoPtr<fvMesh> baseMeshPtr = fvMeshTools::newMesh
            (
                IOobject
                (
                    regionName,
                    baseRunTime.timeName(),
                    baseRunTime,
                    IOobject::MUST_READ
                ),
                true            // read on master only
            );

            setBasicGeometry(baseMeshPtr());


            Info<< "Reading local, decomposed mesh" << endl;
            autoPtr<fvMesh> meshPtr = loadOrCreateMesh
            (
                decompose,
                IOobject
                (
                    regionName,
                    baseMeshPtr().facesInstance(),
                    runTime,
                    Foam::IOobject::MUST_READ
                )
            );
            fvMesh& mesh = meshPtr();

            if (writeHandler.valid() && Pstream::master())
            {
                // Remove any left-over empty processor directories created
                // by loadOrCreateMesh to get around e.g. collated start-up
                // problems. Should not happen in reconstruct mode ...
                const bool oldParRun = Pstream::parRun(false);
                removeEmptyDirs(mesh.time().path());
                Pstream::parRun(oldParRun);
            }


            // Read addressing back to base mesh
            autoPtr<mapDistributePolyMesh> distMap;
            readProcAddressing(mesh, baseMeshPtr, distMap);

            // Construct field mapper
            autoPtr<parFvFieldReconstructor> fvReconstructorPtr
            (
                new parFvFieldReconstructor
                (
                    baseMeshPtr(),
                    mesh,
                    distMap(),
                    Pstream::master()       // do I need to write?
                )
            );



            // Since we start from Times[0] and not runTime.timeName() we
            // might overlook point motion in the first timestep
            // (since mesh.readUpdate() below will not be triggered). Instead
            // detect points by hand
            if (mesh.pointsInstance() != mesh.facesInstance())
            {
                Info<< "    Detected initial mesh motion;"
                    << " reconstructing points" << nl
                    << endl;
                fvReconstructorPtr().reconstructPoints();
            }


            // Loop over all times
            forAll(timeDirs, timeI)
            {
                if (newTimes && masterTimeDirSet.found(timeDirs[timeI].name()))
                {
                    Info<< "Skipping time " << timeDirs[timeI].name()
                        << endl << endl;
                    continue;
                }

                // Set time for global database
                runTime.setTime(timeDirs[timeI], timeI);
                baseRunTime.setTime(timeDirs[timeI], timeI);

                Info<< "Time = " << runTime.timeName() << endl << endl;


                // Check if any new meshes need to be read.
                fvMesh::readUpdateState procStat = mesh.readUpdate();

                if (procStat == fvMesh::POINTS_MOVED)
                {
                    Info<< "    Dected mesh motion; reconstructing points" << nl
                        << endl;
                    fvReconstructorPtr().reconstructPoints();
                }
                else if
                (
                    procStat == fvMesh::TOPO_CHANGE
                 || procStat == fvMesh::TOPO_PATCH_CHANGE
                )
                {
                    Info<< "    Detected topology change;"
                        << " reconstructing addressing" << nl << endl;

                    if (baseMeshPtr)
                    {
                        // Cannot do a baseMesh::readUpdate() since not all
                        // processors will have mesh files. So instead just
                        // recreate baseMesh
                        baseMeshPtr.clear();
                        baseMeshPtr = fvMeshTools::newMesh
                        (
                            IOobject
                            (
                                regionName,
                                baseRunTime.timeName(),
                                baseRunTime,
                                IOobject::MUST_READ
                            ),
                            true            // read on master only
                        );
                    }

                    // Re-read procXXXaddressing
                    readProcAddressing(mesh, baseMeshPtr, distMap);

                    // Reset field mapper
                    fvReconstructorPtr.reset
                    (
                        new parFvFieldReconstructor
                        (
                            baseMeshPtr(),
                            mesh,
                            distMap(),
                            Pstream::master()
                        )
                    );
                    lagrangianReconstructorPtr.clear();
                }


                // Get list of objects
                IOobjectList objects(mesh, runTime.timeName());


                // Mesh fields (vol, surface, volInternal)
                reconstructMeshFields
                (
                    fvReconstructorPtr(),
                    objects,
                    selectedFields
                );

                // Clouds (note: might not be present on all processors)
                reconstructLagrangian
                (
                    lagrangianReconstructorPtr,
                    baseMeshPtr(),
                    mesh,
                    distMap(),
                    selectedLagrangianFields
                );

                // If there are any "uniform" directories copy them from
                // the master processor
                copyUniform
                (
                    fh,
                    decompose,
                    reconstruct,
                    mesh.time().timeName(),
                    mesh,
                    baseMeshPtr()
                );
                // Non-region specific. Note: should do outside region loop
                // but would then have to replicate the whole time loop ...
                copyUniform
                (
                    fh,
                    decompose,
                    reconstruct,
                    mesh.time().timeName(),
                    mesh.time(),                // runTime
                    baseMeshPtr().time()        // baseRunTime
                );
            }
        }
    }
    else
    {
        // decompose or redistribution mode.
        //  decompose : master : read from parent dir
        //              slave  : dummy mesh
        //  redistribute : all read mesh or dummy mesh

        // Time coming from processor0 (or undecomposed if no processor0)
        scalar masterTime;
        if (decompose)
        {
            // Use base time. This is to handle e.g. startTime = latestTime
            // which will not do anything if there are no processor directories
            masterTime = timeSelector::selectIfPresent
            (
                baseRunTime,
                args
            )[0].value();
        }
        else
        {
            masterTime = timeSelector::selectIfPresent
            (
                runTime,
                args
            )[0].value();
        }
        Pstream::scatter(masterTime);
        Info<< "Setting time to that of master or undecomposed case : "
            << masterTime << endl;
        runTime.setTime(masterTime, 0);
        baseRunTime.setTime(masterTime, 0);

        // Save old time name (since might be incremented)
        const word oldTimeName(runTime.timeName());

        forAll(regionNames, regioni)
        {
            const word& regionName = regionNames[regioni];
            const fileName meshSubDir
            (
                regionName == polyMesh::defaultRegion
              ? fileName(polyMesh::meshSubDir)
              : regionNames[regioni]/polyMesh::meshSubDir
            );

            if (decompose)
            {
                Info<< "\n\nDecomposing mesh " << regionName << nl << endl;
            }
            else
            {
                Info<< "\n\nRedistributing mesh " << regionName << nl << endl;
            }


            // Get time instance directory
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // At this point we should be able to read at least a mesh on
            // processor0. Note the changing of the processor0 casename to
            // enforce it to read/write from the undecomposed case

            fileName masterInstDir;
            if (Pstream::master())
            {
                if (decompose)
                {
                    Info<< "Setting caseName to " << baseRunTime.caseName()
                        << " to find undecomposed mesh" << endl;
                    runTime.caseName() = baseRunTime.caseName();
                    runTime.processorCase(false);
                }

                const bool oldParRun = Pstream::parRun(false);
                masterInstDir = runTime.findInstance
                (
                    meshSubDir,
                    "faces",
                    IOobject::READ_IF_PRESENT
                );
                Pstream::parRun(oldParRun);

                if (decompose)
                {
                    Info<< "Restoring caseName to " << proc0CaseName << endl;
                    runTime.caseName() = proc0CaseName;
                    runTime.processorCase(oldProcCase);
                }
            }
            Pstream::scatter(masterInstDir);

            // Check who has a mesh
            const fileName meshPath(runTime.path()/masterInstDir/meshSubDir);
            const boolList haveMesh(haveFacesFile(meshPath));

            // Collect objectPath of polyMesh for the current file handler. This
            // is where the mesh would be written if it didn't exist already.
            fileNameList meshDir(Pstream::nProcs());
            {
                const fileName fName
                (
                    fileHandler().objectPath
                    (
                        IOobject("faces", masterInstDir/meshSubDir, runTime),
                        word::null
                    )
                );
                meshDir[Pstream::myProcNo()] = fName.path();
                Pstream::gatherList(meshDir);
                Pstream::scatterList(meshDir);
                //Info<< "Per processor faces dirs:" << nl
                //    << "    " << meshDir << nl << endl;
            }


            // Load mesh (or create dummy one)
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            if (Pstream::master() && decompose)
            {
                Info<< "Setting caseName to " << baseRunTime.caseName()
                    << " to read undecomposed mesh" << endl;
                runTime.caseName() = baseRunTime.caseName();
                runTime.processorCase(false);
            }

            autoPtr<fvMesh> meshPtr = loadOrCreateMesh
            (
                decompose,
                //haveMesh[Pstream::myProcNo()],
                IOobject
                (
                    regionName,
                    masterInstDir,
                    runTime,
                    Foam::IOobject::MUST_READ
                )
            );
            fvMesh& mesh = meshPtr();

            if (writeHandler.valid())
            {
                // Remove any left-over empty processor directories created
                // by loadOrCreateMesh to get around the collated start-up
                // problems
                if (Pstream::master())  //fileHandler().comm()))
                {
                    const auto myProci = UPstream::myProcNo();  //comm()
                    const auto& procs = UPstream::procID
                    (
                        UPstream::worldComm
                    );
                    const bool oldParRun = Pstream::parRun(false);
                    for (const auto proci : procs)
                    {
                        if
                        (
                           !haveMesh[proci]
                         && meshDir[proci] != meshDir[myProci]
                        )
                        {
                            Info<< "Deleting mesh dir:" << meshDir[proci]
                                << endl;
                            rmDir(meshDir[proci]);
                        }
                    }

                    // Remove empty directory
                    removeEmptyDirs(mesh.time().path());

                    Pstream::parRun(oldParRun);
                }
            }


            if (Pstream::master() && decompose)
            {
                Info<< "Restoring caseName to " << proc0CaseName << endl;
                runTime.caseName() = proc0CaseName;
                runTime.processorCase(oldProcCase);
            }

            const label nOldCells = mesh.nCells();
            //Pout<< "Loaded mesh : nCells:" << nOldCells
            //    << " nPatches:" << mesh.boundaryMesh().size() << endl;


            // Determine decomposition
            // ~~~~~~~~~~~~~~~~~~~~~~~

            label nDestProcs;
            labelList finalDecomp;
            determineDecomposition
            (
                baseRunTime,
                decompDictFile,
                decompose,
                proc0CaseName,
                mesh,
                writeCellDist,

                nDestProcs,
                finalDecomp
            );

            if (dryrun)
            {
                if (!Pstream::master() && !haveMesh[Pstream::myProcNo()])
                {
                    // Remove dummy mesh created by loadOrCreateMesh
                    const bool oldParRun = Pstream::parRun(false);
                    mesh.removeFiles();
                    rmDir(mesh.objectRegistry::objectPath());
                    Pstream::parRun(oldParRun);  // Restore parallel state
                }
                continue;
            }



            wordList cloudNames;
            List<wordList> fieldNames;

            // Detect lagrangian fields
            if (Pstream::master() && decompose)
            {
                runTime.caseName() = baseRunTime.caseName();
                runTime.processorCase(false);
            }
            parLagrangianRedistributor::findClouds
            (
                mesh,
                cloudNames,
                fieldNames
            );

            // Read lagrangian fields and store on cloud (objectRegistry)
            PtrList<unmappedPassivePositionParticleCloud> clouds
            (
                cloudNames.size()
            );
            readLagrangian
            (
                mesh,
                cloudNames,
                selectedLagrangianFields,
                clouds
            );
            if (Pstream::master() && decompose)
            {
                runTime.caseName() = proc0CaseName;
                runTime.processorCase(oldProcCase);
            }

            // Load fields, do all distribution (mesh and fields - but not
            // lagrangian fields; these are done later)
            autoPtr<mapDistributePolyMesh> distMap = redistributeAndWrite
            (
                std::move(writeHandler),
                baseRunTime,
                haveMesh,
                meshSubDir,
                true,           // read fields
                decompose,      // decompose, i.e. read from undecomposed case
                false,          // no reconstruction
                overwrite,
                proc0CaseName,
                nDestProcs,
                finalDecomp,
                masterInstDir,
                mesh
            );


            // Redistribute any clouds
            redistributeLagrangian
            (
                lagrangianReconstructorPtr,
                mesh,
                nOldCells,
                distMap(),
                clouds
            );


            // Copy region-specific uniform
            // (e.g. solid/uniform/cumulativeContErr)
            copyUniform
            (
                fh,
                decompose,
                reconstruct,
                oldTimeName,    // provided read time
                mesh,           // read location is mesh (but oldTimeName)
                mesh            // write location is mesh
            );
        }

        // Copy non-region specific uniform (e.g. uniform/time)
        copyUniform
        (
            fh,
            decompose,
            reconstruct,
            oldTimeName,    // provided read time
            (               // read location
                decompose
              ? baseRunTime
              : runTime
            ),
            runTime         // writing location
        );
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
