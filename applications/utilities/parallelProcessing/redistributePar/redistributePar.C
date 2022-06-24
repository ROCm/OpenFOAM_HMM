/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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
#include "fieldsDistributor.H"
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

#include "parFvFieldDistributor.H"
#include "parPointFieldDistributor.H"
#include "hexRef8Data.H"
#include "meshRefinement.H"
#include "pointFields.H"

#include "faMeshSubset.H"
#include "faMeshTools.H"
#include "faMeshDistributor.H"
#include "parFaFieldDistributorCache.H"

#include "redistributeLagrangian.H"

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
            Foam::mkDir(timePath);
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


void printMeshData(const polyMesh& mesh)
{
    // Collect all data on master

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

    const globalIndex globalCells(mesh.nCells());
    const globalIndex globalBoundaryFaces(mesh.nBoundaryFaces());

    label maxProcCells = 0;
    label maxProcFaces = 0;
    label totProcFaces = 0;
    label maxProcPatches = 0;
    label totProcPatches = 0;

    for (const int proci : Pstream::allProcs())
    {
        const label nLocalCells = globalCells.localSize(proci);
        const label nBndFaces = globalBoundaryFaces.localSize(proci);

        Info<< nl
            << "Processor " << proci;

        if (!nLocalCells)
        {
            Info<< " (empty)" << endl;
            continue;
        }
        else
        {
            Info<< nl
                << "    Number of cells = " << nLocalCells << endl;
        }

        label nProcFaces = 0;
        const labelList& nei = patchNeiProcNo[proci];

        forAll(patchNeiProcNo[proci], i)
        {
            Info<< "    Number of faces shared with processor "
                << patchNeiProcNo[proci][i] << " = "
                << patchSize[proci][i] << nl;

            nProcFaces += patchSize[proci][i];
        }

        {
            Info<< "    Number of processor patches = " << nei.size() << nl
                << "    Number of processor faces = " << nProcFaces << nl
                << "    Number of boundary faces = "
                << nBndFaces-nProcFaces << endl;
        }

        maxProcCells = max(maxProcCells, nLocalCells);
        totProcFaces += nProcFaces;
        totProcPatches += nei.size();
        maxProcFaces = max(maxProcFaces, nProcFaces);
        maxProcPatches = max(maxProcPatches, nei.size());
    }

    // Summary stats

    Info<< nl
        << "Number of processor faces = " << (totProcFaces/2) << nl
        << "Max number of cells = " << maxProcCells;

    if (maxProcCells != globalCells.totalSize())
    {
        scalar avgValue = scalar(globalCells.totalSize())/Pstream::nProcs();

        Info<< " (" << 100.0*(maxProcCells-avgValue)/avgValue
            << "% above average " << avgValue << ')';
    }
    Info<< nl;

    Info<< "Max number of processor patches = " << maxProcPatches;
    if (totProcPatches)
    {
        scalar avgValue = scalar(totProcPatches)/Pstream::nProcs();

        Info<< " (" << 100.0*(maxProcPatches-avgValue)/avgValue
            << "% above average " << avgValue << ')';
    }
    Info<< nl;

    Info<< "Max number of faces between processors = " << maxProcFaces;
    if (totProcFaces)
    {
        scalar avgValue = scalar(totProcFaces)/Pstream::nProcs();

        Info<< " (" << 100.0*(maxProcFaces-avgValue)/avgValue
            << "% above average " << avgValue << ')';
    }
    Info<< nl << endl;
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
            IOobject::NO_WRITE,
            false                   // do not register
        ),
        mesh,
        dimensionedScalar(name, dimless, -1),
        zeroGradientFvPatchScalarField::typeName
    );

    forAll(procCells, celli)
    {
        procCells[celli] = decomp[celli];
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
            << "You have selected decomposition method \""
            << decomposer.type() << "\n"
            << "    which does not synchronise decomposition across"
               " processor patches.\n"
               "    You might want to select a decomposition method"
               " that is aware of this. Continuing...." << endl;
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
        Info<< "Restoring caseName" << endl;
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
                Info<< "Restoring caseName" << endl;
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


// Variant of GeometricField::correctBoundaryConditions that only
// evaluates selected patch fields
template<class CoupledPatchType, class GeoField>
void correctCoupledBoundaryConditions(fvMesh& mesh)
{
    for (GeoField& fld : mesh.sorted<GeoField>())
    {
        fld.boundaryFieldRef().template evaluateCoupled<CoupledPatchType>();
    }
}


// Inplace redistribute mesh and any fields
autoPtr<mapDistributePolyMesh> redistributeAndWrite
(
    autoPtr<fileOperation>&& writeHandler,
    const Time& baseRunTime,
    const fileName& proc0CaseName,

    // Controls
    const bool doReadFields,
    const bool decompose,       // decompose, i.e. read from undecomposed case
    const bool reconstruct,
    const bool overwrite,

    // Decomposition information
    const label nDestProcs,
    const labelList& decomp,

    // Mesh information
    const boolList& volMeshOnProc,
    const fileName& volMeshInstance,
    fvMesh& mesh
)
{
    Time& runTime = const_cast<Time&>(mesh.time());
    const bool oldProcCase = runTime.processorCase();

    //// Print some statistics
    //Info<< "Before distribution:" << endl;
    //printMeshData(mesh);

    // Storage of fields

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

    PtrList<pointScalarField> pointScalarFields;
    PtrList<pointVectorField> pointVectorFields;
    PtrList<pointTensorField> pointTensorFields;
    PtrList<pointSphericalTensorField> pointSphTensorFields;
    PtrList<pointSymmTensorField> pointSymmTensorFields;

    // Self-contained pointMesh for reading pointFields
    const pointMesh oldPointMesh(mesh);

    // Track how many (if any) pointFields are read/mapped
    label nPointFields = 0;

    parPointFieldDistributor pointDistributor
    (
        oldPointMesh,   // source mesh
        false,          // savePoints=false (ie, delay until later)
        false           // Do not write
    );


    if (doReadFields)
    {
        // Create 0 sized mesh to do all the generation of zero sized
        // fields on processors that have zero sized meshes. Note that this is
        // only necessary on master but since polyMesh construction with
        // Pstream::parRun does parallel comms we have to do it on all
        // processors
        autoPtr<fvMeshSubset> subsetterPtr;

        // Missing a volume mesh somewhere?
        if (volMeshOnProc.found(false))
        {
            // A zero-sized mesh with boundaries.
            // This is used to create zero-sized fields.
            subsetterPtr.reset(new fvMeshSubset(mesh, zero{}));
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
            << " mesh:" << mesh.objectRegistry::objectRelPath()
            << " have objects:" << objects.names() << endl;

        // We don't want to map the decomposition (mapping already tested when
        // mapping the cell centre field)
        auto iter = objects.find("cellDist");
        if (iter.found())
        {
            objects.erase(iter);
        }


        if (Pstream::master() && decompose)
        {
            runTime.caseName() = baseRunTime.caseName();
            runTime.processorCase(false);
        }

        // Field reading

        #undef  doFieldReading
        #define doFieldReading(Storage)                                       \
        {                                                                     \
            fieldsDistributor::readFields                                     \
            (                                                                 \
                volMeshOnProc, mesh, subsetterPtr, objects, Storage           \
            );                                                                \
        }

        // volField
        doFieldReading(volScalarFields);
        doFieldReading(volVectorFields);
        doFieldReading(volSphereTensorFields);
        doFieldReading(volSymmTensorFields);
        doFieldReading(volTensorFields);

        // surfaceField
        doFieldReading(surfScalarFields);
        doFieldReading(surfVectorFields);
        doFieldReading(surfSphereTensorFields);
        doFieldReading(surfSymmTensorFields);
        doFieldReading(surfTensorFields);

        // Dimensioned internal fields
        doFieldReading(dimScalarFields);
        doFieldReading(dimVectorFields);
        doFieldReading(dimSphereTensorFields);
        doFieldReading(dimSymmTensorFields);
        doFieldReading(dimTensorFields);

        // pointFields
        nPointFields = 0;

        #undef  doFieldReading
        #define doFieldReading(Storage)                                       \
        {                                                                     \
            fieldsDistributor::readFields                                     \
            (                                                                 \
                volMeshOnProc, oldPointMesh, subsetterPtr, objects, Storage,  \
                true  /* (deregister field) */                                \
            );                                                                \
            nPointFields += Storage.size();                                   \
        }

        doFieldReading(pointScalarFields);
        doFieldReading(pointVectorFields);
        doFieldReading(pointSphTensorFields);
        doFieldReading(pointSymmTensorFields);
        doFieldReading(pointTensorFields);
        #undef doFieldReading


        // Done reading

        if (Pstream::master() && decompose)
        {
            runTime.caseName() = proc0CaseName;
            runTime.processorCase(oldProcCase);
        }
    }

    // Save pointMesh information before any topology changes occur!
    if (nPointFields)
    {
        pointDistributor.saveMeshPoints();
    }


    // Mesh distribution engine
    fvMeshDistribute distributor(mesh);

    // Do all the distribution of mesh and fields
    autoPtr<mapDistributePolyMesh> distMap = distributor.distribute(decomp);

    // Print some statistics
    Info<< "After distribution:" << endl;
    printMeshData(mesh);

    // Get other side of processor boundaries
    do
    {
        #undef  doCorrectCoupled
        #define doCorrectCoupled(FieldType)  \
        correctCoupledBoundaryConditions<processorFvPatch, FieldType>(mesh);

        doCorrectCoupled(volScalarField);
        doCorrectCoupled(volVectorField);
        doCorrectCoupled(volSphericalTensorField);
        doCorrectCoupled(volSymmTensorField);
        doCorrectCoupled(volTensorField);
        #undef doCorrectCoupled
    }
    while (false);

    // No update surface fields


    // Map pointFields
    if (nPointFields)
    {
        // Construct new pointMesh from distributed mesh
        const pointMesh& newPointMesh = pointMesh::New(mesh);

        pointDistributor.resetTarget(newPointMesh, distMap());

        pointDistributor.distributeAndStore(pointScalarFields);
        pointDistributor.distributeAndStore(pointVectorFields);
        pointDistributor.distributeAndStore(pointSphTensorFields);
        pointDistributor.distributeAndStore(pointSymmTensorFields);
        pointDistributor.distributeAndStore(pointTensorFields);
    }


    // Set the minimum write precision
    IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));


    if (!overwrite)
    {
        ++runTime;
        mesh.setInstance(runTime.timeName());
    }
    else
    {
        mesh.setInstance(volMeshInstance);
    }


    // Register mapDistributePolyMesh for automatic writing...
    IOmapDistributePolyMeshRef distMapRef
    (
        IOobject
        (
            "procAddressing",
            mesh.facesInstance(),
            polyMesh::meshSubDir,
            mesh.thisDb(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        distMap()
    );


    if (reconstruct)
    {
        if (Pstream::master())
        {
            Info<< "Setting caseName to " << baseRunTime.caseName()
                << " to write reconstructed mesh (and fields)." << endl;
            runTime.caseName() = baseRunTime.caseName();
            const bool oldProcCase(runTime.processorCase(false));

            mesh.write();
            topoSet::removeFiles(mesh);

            // Now we've written all. Reset caseName on master
            Info<< "Restoring caseName" << endl;
            runTime.caseName() = proc0CaseName;
            runTime.processorCase(oldProcCase);
        }
    }
    else
    {
        autoPtr<fileOperation> defaultHandler;
        if (writeHandler)
        {
            defaultHandler = fileHandler(std::move(writeHandler));
        }

        mesh.write();

        if (defaultHandler)
        {
            writeHandler = fileHandler(std::move(defaultHandler));
        }
        topoSet::removeFiles(mesh);
    }
    Info<< "Written redistributed mesh to "
        << mesh.facesInstance() << nl << endl;


    if (decompose || reconstruct)
    {
        // Decompose (1 -> N) or reconstruct (N -> 1)
        // so {boundary,cell,face,point}ProcAddressing have meaning
        fvMeshTools::writeProcAddressing
        (
            mesh,
            distMap(),
            decompose,
            std::move(writeHandler)
        );
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
        refData.distribute(distMap());


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
                Info<< "Restoring caseName" << endl;
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
    //        cellSets[i].distribute(distMap());
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
    //                cellSets[i].distribute(distMap());
    //            }
    //
    //            // Now we've written all. Reset caseName on master
    //            Info<< "Restoring caseName" << endl;
    //            runTime.caseName() = proc0CaseName;
    //            runTime.processorCase(oldProcCase);
    //        }
    //    }
    //    else
    //    {
    //        forAll(cellSets, i)
    //        {
    //            cellSets[i].distribute(distMap());
    //        }
    //    }
    //}


    return distMap;
}


/*---------------------------------------------------------------------------*\
                                     main
\*---------------------------------------------------------------------------*/

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
    argList::addVerboseOption("Additional verbosity");
    argList::addDryRunOption
    (
        "Test without writing the decomposition. "
        "Changes -cellDist to only write volScalarField."
    );
    argList::addVerboseOption("Additional verbosity");
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
    argList::addVerboseOption
    (
        "Additional verbosity. (Can be used multiple times)"
    );
    argList::addBoolOption
    (
        "no-finite-area",
        "Suppress finiteArea mesh/field handling",
        true  // Advanced option
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
        writeHandler
      ? writeHandler()
      : fileHandler()
    );

    // Make sure to call findTimes on all processors to force caching of
    // time directories
    (void)fileHandler().findTimes(args.path(), "constant");

    // Need this line since we don't include "setRootCase.H"
    #include "foamDlOpenLibs.H"

    const bool reconstruct = args.found("reconstruct");
    const bool writeCellDist = args.found("cellDist");
    const bool dryrun = args.dryRun();
    const bool newTimes = args.found("newTimes");
    const int optVerbose = args.verbose();

    const bool doFiniteArea = !args.found("no-finite-area");
    bool decompose = args.found("decompose");
    bool overwrite = args.found("overwrite");

    if (optVerbose)
    {
        // Report on output
        faMeshDistributor::verbose_ = 1;
        parPointFieldDistributor::verbose_ = 1;
    }

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

    if ((decompose || reconstruct) && !overwrite)
    {
        overwrite = true;
        WarningInFunction
            << "Working in -decompose or -reconstruct mode:"
               " automatically implies -overwrite" << nl << endl;
    }

    if (!Pstream::parRun())
    {
        FatalErrorInFunction
            << ": This utility can only be run parallel"
            << exit(FatalError);
    }


    if (!Foam::isDir(args.rootPath()))
    {
        FatalErrorInFunction
            << ": cannot open root directory " << args.rootPath()
            << exit(FatalError);
    }

    // Detect if running data-distributed (multiple roots)
    bool nfs = true;
    {
        List<fileName> roots(1, args.rootPath());
        Pstream::combineAllGather(roots, ListOps::uniqueEqOp<fileName>());
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
    if (Foam::isDir(procDir))
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


    Info<< "Create undecomposed database" << nl << endl;
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
    autoPtr<parLagrangianDistributor> lagrangianDistributorPtr;


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
            << "Reconstructing mesh and addressing" << nl << endl;

        forAll(regionNames, regioni)
        {
            const word& regionName = regionNames[regioni];
            const word& regionDir = polyMesh::regionName(regionName);

            const fileName volMeshSubDir(regionDir/polyMesh::meshSubDir);
            const fileName areaMeshSubDir(regionDir/faMesh::meshSubDir);

            Info<< nl
                << "Reconstructing mesh " << regionDir << nl << endl;

            bool areaMeshDetected = false;

            // Loop over all times
            forAll(timeDirs, timeI)
            {
                // Set time for global database
                runTime.setTime(timeDirs[timeI], timeI);
                baseRunTime.setTime(timeDirs[timeI], timeI);

                Info<< "Time = " << runTime.timeName() << endl << endl;

                // Where meshes are
                fileName volMeshInstance;
                fileName areaMeshInstance;

                volMeshInstance = runTime.findInstance
                (
                    volMeshSubDir,
                    "faces",
                    IOobject::READ_IF_PRESENT
                );

                if (doFiniteArea)
                {
                    areaMeshInstance = runTime.findInstance
                    (
                        areaMeshSubDir,
                        "faceLabels",
                        IOobject::READ_IF_PRESENT
                    );
                }

                Pstream::broadcasts
                (
                    UPstream::worldComm,
                    volMeshInstance,
                    areaMeshInstance
                );


                // Check processors have meshes
                // - check for 'faces' file (polyMesh)
                // - check for 'faceLabels' file (faMesh)
                boolList volMeshOnProc;
                boolList areaMeshOnProc;

                volMeshOnProc = haveMeshFile
                (
                    runTime,
                    volMeshInstance/volMeshSubDir,
                    "faces"
                );

                if (doFiniteArea)
                {
                    areaMeshOnProc = haveMeshFile
                    (
                        runTime,
                        areaMeshInstance/areaMeshSubDir,
                        "faceLabels"
                    );

                    areaMeshDetected = areaMeshOnProc.found(true);
                }


                // Addressing back to reconstructed mesh as xxxProcAddressing.
                // - all processors have consistent faceProcAddressing
                // - processors without a mesh don't need faceProcAddressing


                // Note: filePath searches up on processors that don't have
                //       processor if instance = constant so explicitly check
                //       found filename.
                bool haveVolAddressing = false;
                if (volMeshOnProc[Pstream::myProcNo()])
                {
                    // Read faces (just to know their size)
                    faceCompactIOList faces
                    (
                        IOobject
                        (
                            "faces",
                            volMeshInstance,
                            volMeshSubDir,
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
                            volMeshInstance,
                            volMeshSubDir,
                            runTime,
                            IOobject::READ_IF_PRESENT
                        )
                    );

                    haveVolAddressing =
                    (
                        faceProcAddressing.headerOk()
                     && faceProcAddressing.size() == faces.size()
                    );
                }
                else
                {
                    // Have no mesh. Don't need addressing
                    haveVolAddressing = true;
                }

                bool haveAreaAddressing = false;
                if (areaMeshOnProc[Pstream::myProcNo()])
                {
                    // Read faces (just to know their size)
                    labelIOList faceLabels
                    (
                        IOobject
                        (
                            "faceLabels",
                            areaMeshInstance,
                            areaMeshSubDir,
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
                            areaMeshInstance,
                            areaMeshSubDir,
                            runTime,
                            IOobject::READ_IF_PRESENT
                        )
                    );

                    haveAreaAddressing =
                    (
                        faceProcAddressing.headerOk()
                     && faceProcAddressing.size() == faceLabels.size()
                    );
                }
                else if (areaMeshDetected)
                {
                    // Have no mesh. Don't need addressing
                    haveAreaAddressing = true;
                }


                // Additionally check for master faces being readable. Could
                // do even more checks, e.g. global number of cells same
                // as cellProcAddressing

                bool volMeshHaveUndecomposed = false;
                bool areaMeshHaveUndecomposed = false;

                if (Pstream::master())
                {
                    Info<< "Checking " << baseRunTime.caseName()
                        << " for undecomposed volume and area meshes..."
                        << endl;

                    const bool oldParRun = Pstream::parRun(false);

                    // Volume
                    {
                        faceCompactIOList facesIO
                        (
                            IOobject
                            (
                                "faces",
                                volMeshInstance,
                                volMeshSubDir,
                                baseRunTime,
                                IOobject::NO_READ
                            ),
                            label(0)
                        );
                        volMeshHaveUndecomposed = facesIO.headerOk();
                    }

                    // Area
                    if (doFiniteArea)
                    {
                        labelIOList labelsIO
                        (
                            IOobject
                            (
                                "faceLabels",
                                areaMeshInstance,
                                areaMeshSubDir,
                                baseRunTime,
                                IOobject::NO_READ
                            )
                        );
                        areaMeshHaveUndecomposed = labelsIO.headerOk();
                    }

                    Pstream::parRun(oldParRun);  // Restore parallel state
                }

                Pstream::broadcasts
                (
                    UPstream::worldComm,
                    volMeshHaveUndecomposed,
                    areaMeshHaveUndecomposed
                );

                // Report
                {
                    Info<< "    volume mesh ["
                        << volMeshHaveUndecomposed << "] : "
                        << volMeshInstance << nl
                        << "    area   mesh ["
                        << areaMeshHaveUndecomposed << "] : "
                        << areaMeshInstance << nl
                        << endl;
                }


                if
                (
                    !volMeshHaveUndecomposed
                 || !returnReduce(haveVolAddressing, andOp<bool>())
                )
                {
                    Info<< "No undecomposed mesh. Creating from: "
                        << volMeshInstance << endl;

                    if (areaMeshHaveUndecomposed)
                    {
                        areaMeshHaveUndecomposed = false;
                        Info<< "Also ignore any undecomposed area mesh"
                            << endl;
                    }

                    autoPtr<fvMesh> volMeshPtr = fvMeshTools::loadOrCreateMesh
                    (
                        IOobject
                        (
                            regionName,
                            volMeshInstance,
                            runTime,
                            Foam::IOobject::MUST_READ
                        ),
                        decompose
                    );
                    fvMeshTools::setBasicGeometry(volMeshPtr());
                    fvMesh& mesh = volMeshPtr();


                    Info<< nl << "Reconstructing mesh" << nl << endl;

                    // Reconstruct (1 processor)
                    const label nDestProcs(1);
                    const labelList finalDecomp(mesh.nCells(), Zero);

                    redistributeAndWrite
                    (
                        std::move(writeHandler),
                        baseRunTime,
                        proc0CaseName,

                        // Controls
                        false,      // do not read fields
                        false,      // do not read undecomposed case on proc0
                        true,       // write redistributed files to proc0
                        overwrite,

                        // Decomposition information
                        nDestProcs,
                        finalDecomp,

                        // For finite-volume
                        volMeshOnProc,
                        volMeshInstance,
                        mesh
                    );
                }


                // Similarly for finiteArea
                // - may or may not have undecomposed mesh
                // - may or may not have decomposed meshes

                if
                (
                    areaMeshOnProc.found(true)  // ie, areaMeshDetected
                 &&
                    (
                        !areaMeshHaveUndecomposed
                     || !returnReduce(haveAreaAddressing, andOp<bool>())
                    )
                )
                {
                    Info<< "Loading area mesh from "
                        << areaMeshInstance << endl;

                    Info<< "    getting volume mesh support" << endl;

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
                    fvMeshTools::setBasicGeometry(baseMeshPtr());

                    autoPtr<fvMesh> volMeshPtr = fvMeshTools::loadOrCreateMesh
                    (
                        IOobject
                        (
                            regionName,
                            baseMeshPtr().facesInstance(),
                            runTime,
                            Foam::IOobject::MUST_READ
                        ),
                        decompose
                    );
                    fvMeshTools::setBasicGeometry(volMeshPtr());
                    fvMesh& mesh = volMeshPtr();

                    // Read volume proc addressing back to base mesh
                    autoPtr<mapDistributePolyMesh> distMap
                    (
                        fvMeshTools::readProcAddressing(mesh, baseMeshPtr)
                    );


                    autoPtr<faMesh> areaMeshPtr = faMeshTools::loadOrCreateMesh
                    (
                        IOobject
                        (
                            regionName,
                            areaMeshInstance,
                            runTime,
                            Foam::IOobject::MUST_READ
                        ),
                        mesh,  // <- The referenced polyMesh (from above)
                        decompose
                    );
                    faMesh& areaMesh = areaMeshPtr();

                    faMeshTools::forceDemandDriven(areaMesh);
                    faMeshTools::unregisterMesh(areaMesh);

                    autoPtr<faMesh> areaBaseMeshPtr;

                    // Reconstruct using polyMesh distribute map
                    mapDistributePolyMesh faDistMap
                    (
                        faMeshDistributor::distribute
                        (
                            areaMesh,
                            distMap(),      // The polyMesh distMap
                            baseMeshPtr(),  // Target polyMesh
                            areaBaseMeshPtr
                        )
                    );

                    faMeshTools::forceDemandDriven(areaBaseMeshPtr());
                    faMeshTools::unregisterMesh(areaBaseMeshPtr());


                    if (Pstream::master())
                    {
                        Info<< "Setting caseName to " << baseRunTime.caseName()
                            << " to write reconstructed area mesh." << endl;
                        runTime.caseName() = baseRunTime.caseName();
                        const bool oldProcCase(runTime.processorCase(false));

                        areaBaseMeshPtr().write();

                        // Now we've written all. Reset caseName on master
                        Info<< "Restoring caseName" << endl;
                        runTime.caseName() = proc0CaseName;
                        runTime.processorCase(oldProcCase);
                    }

                    // Update for the reconstructed procAddressing
                    faMeshTools::writeProcAddressing
                    (
                        areaBaseMeshPtr(),  // Reconstruct location
                        faDistMap,
                        false,              // decompose=false
                        std::move(writeHandler),
                        areaMeshPtr.get()   // procMesh
                    );
                }
            }

            // Make sure all is finished writing until re-reading in pass2
            // below
            fileHandler().flush();


            // Pass2 : read mesh and addressing and reconstruct fields
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            Info<< nl
                << "Reconstructing fields" << nl << endl;

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

            fvMeshTools::setBasicGeometry(baseMeshPtr());

            Info<< "Reading local, decomposed mesh" << endl;
            autoPtr<fvMesh> volMeshPtr = fvMeshTools::loadOrCreateMesh
            (
                IOobject
                (
                    regionName,
                    baseMeshPtr().facesInstance(),
                    runTime,
                    Foam::IOobject::MUST_READ
                ),
                decompose
            );
            fvMesh& mesh = volMeshPtr();


            // Similarly for finiteArea
            autoPtr<faMesh> areaBaseMeshPtr;
            autoPtr<faMesh> areaMeshPtr;
            autoPtr<faMeshDistributor> faDistributor;
            mapDistributePolyMesh areaDistMap;

            if (areaMeshDetected)
            {
                areaBaseMeshPtr = faMeshTools::newMesh
                (
                    IOobject
                    (
                        regionName,
                        baseRunTime.timeName(),
                        baseRunTime,
                        IOobject::MUST_READ
                    ),
                    baseMeshPtr(),
                    true            // read on master only
                );

                areaMeshPtr = faMeshTools::loadOrCreateMesh
                (
                    IOobject
                    (
                        regionName,
                        areaBaseMeshPtr().facesInstance(),
                        runTime,
                        IOobject::MUST_READ
                    ),
                    mesh,
                    decompose
                );

                areaDistMap =
                    faMeshTools::readProcAddressing
                    (
                        areaMeshPtr(),
                        areaBaseMeshPtr
                    );

                faMeshTools::forceDemandDriven(areaMeshPtr());

                // Create an appropriate field distributor
                faDistributor.reset
                (
                    new faMeshDistributor
                    (
                        areaMeshPtr(),      // source
                        areaBaseMeshPtr(),  // target
                        areaDistMap,
                        Pstream::master()   // only write on master
                    )
                );
                // Report some messages. Tbd.
                faMeshDistributor::verbose_ = 1;
            }


            if (writeHandler && Pstream::master())
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
            distMap = fvMeshTools::readProcAddressing(mesh, baseMeshPtr);

            // Construct field mapper
            auto fvDistributorPtr =
                autoPtr<parFvFieldDistributor>::New
                (
                    mesh,            // source
                    baseMeshPtr(),   // target
                    distMap(),
                    UPstream::master()  // Write reconstructed on master
                );

            // Construct point field mapper
            const auto& basePointMesh = pointMesh::New(baseMeshPtr());
            const auto& procPointMesh = pointMesh::New(mesh);

            auto pointFieldDistributorPtr =
                autoPtr<parPointFieldDistributor>::New
                (
                    procPointMesh,   // source
                    basePointMesh,   // target
                    distMap(),
                    false,           // delay
                    UPstream::master()  // Write reconstructed on master
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
                fvDistributorPtr().reconstructPoints();
            }


            // Loop over all times
            forAll(timeDirs, timeI)
            {
                if (newTimes && masterTimeDirSet.found(timeDirs[timeI].name()))
                {
                    Info<< "Skipping time " << timeDirs[timeI].name()
                        << nl << endl;
                    continue;
                }

                // Set time for global database
                runTime.setTime(timeDirs[timeI], timeI);
                baseRunTime.setTime(timeDirs[timeI], timeI);

                Info<< "Time = " << runTime.timeName() << endl << endl;


                // Check if any new meshes need to be read.
                polyMesh::readUpdateState procStat = mesh.readUpdate();

                if (procStat == polyMesh::POINTS_MOVED)
                {
                    Info<< "    Detected mesh motion; reconstructing points"
                        << nl << endl;
                    fvDistributorPtr().reconstructPoints();
                }
                else if
                (
                    procStat == polyMesh::TOPO_CHANGE
                 || procStat == polyMesh::TOPO_PATCH_CHANGE
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

                    // Re-read procAddressing
                    distMap =
                        fvMeshTools::readProcAddressing(mesh, baseMeshPtr);

                    // Reset field mappers

                    fvDistributorPtr.reset
                    (
                        new parFvFieldDistributor
                        (
                            mesh,           // source
                            baseMeshPtr(),  // target
                            distMap(),
                            UPstream::master()  // Write reconstruct on master
                        )
                    );

                    // Construct point field mapper
                    const auto& basePointMesh = pointMesh::New(baseMeshPtr());
                    const auto& procPointMesh = pointMesh::New(mesh);

                    pointFieldDistributorPtr.reset
                    (
                        new parPointFieldDistributor
                        (
                            procPointMesh,  // source
                            basePointMesh,  // target
                            distMap(),
                            false,          // delay until later
                            UPstream::master()  // Write reconstruct on master
                        )
                    );

                    lagrangianDistributorPtr.reset();

                    if (areaMeshPtr)
                    {
                        Info<< "    Discarding finite-area addressing"
                            << " (TODO)" << nl << endl;

                        areaBaseMeshPtr.reset();
                        areaMeshPtr.reset();
                        faDistributor.reset();
                        areaDistMap.clear();
                    }
                }


                // Get list of objects
                IOobjectList objects(mesh, runTime.timeName());


                // Mesh fields (vol, surface, volInternal)
                fvDistributorPtr()
                    .distributeAllFields(objects, selectedFields);

                // pointfields
                // - distribute and write (verbose)
                pointFieldDistributorPtr()
                    .distributeAllFields(objects, selectedFields);


                // Clouds (note: might not be present on all processors)
                reconstructLagrangian
                (
                    lagrangianDistributorPtr,
                    baseMeshPtr(),
                    mesh,
                    distMap(),
                    selectedLagrangianFields
                );

                if (faDistributor)
                {
                    faDistributor()
                        .distributeAllFields(objects, selectedFields);
                }

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
            const word& regionDir = polyMesh::regionName(regionName);

            const fileName volMeshSubDir(regionDir/polyMesh::meshSubDir);
            const fileName areaMeshSubDir(regionDir/faMesh::meshSubDir);

            Info<< nl << nl
                << (decompose ? "Decomposing" : "Redistributing")
                << " mesh " << regionDir << nl << endl;


            // Get time instance directory
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // At this point we should be able to read at least a mesh on
            // processor0. Note the changing of the processor0 casename to
            // enforce it to read/write from the undecomposed case

            fileName volMeshMasterInstance;
            fileName areaMeshMasterInstance;

            // Assume to be true
            bool volMeshHaveUndecomposed = true;
            bool areaMeshHaveUndecomposed = doFiniteArea;

            if (Pstream::master())
            {
                if (decompose)
                {
                    Info<< "Checking undecomposed mesh in case: "
                        << baseRunTime.caseName() << endl;
                    runTime.caseName() = baseRunTime.caseName();
                    runTime.processorCase(false);
                }

                const bool oldParRun = Pstream::parRun(false);
                volMeshMasterInstance = runTime.findInstance
                (
                    volMeshSubDir,
                    "faces",
                    IOobject::READ_IF_PRESENT
                );

                if (doFiniteArea)
                {
                    areaMeshMasterInstance = runTime.findInstance
                    (
                        areaMeshSubDir,
                        "faceLabels",
                        IOobject::READ_IF_PRESENT
                    );

                    // Note: findInstance returns "constant" even if not found,
                    // so recheck now for a false positive.

                    if ("constant" == areaMeshMasterInstance)
                    {
                        const boolList areaMeshOnProc
                        (
                            haveMeshFile
                            (
                                runTime,
                                areaMeshMasterInstance/areaMeshSubDir,
                                "faceLabels",
                                false  // verbose=false
                            )
                        );

                        if (areaMeshOnProc.empty() || !areaMeshOnProc[0])
                        {
                            areaMeshHaveUndecomposed = false;
                        }
                    }
                }

                Pstream::parRun(oldParRun);  // Restore parallel state

                if (decompose)
                {
                    Info<< "    volume mesh ["
                        << volMeshHaveUndecomposed << "] : "
                        << volMeshMasterInstance << nl
                        << "    area   mesh ["
                        << areaMeshHaveUndecomposed << "] : "
                        << areaMeshMasterInstance << nl
                        << nl << nl;

                    // Restoring caseName
                    runTime.caseName() = proc0CaseName;
                    runTime.processorCase(oldProcCase);
                }
            }

            Pstream::broadcasts
            (
                UPstream::worldComm,
                volMeshHaveUndecomposed,
                areaMeshHaveUndecomposed,
                volMeshMasterInstance,
                areaMeshMasterInstance
            );

            // Check processors have meshes
            // - check for 'faces' file (polyMesh)
            // - check for 'faceLabels' file (faMesh)
            boolList volMeshOnProc;
            boolList areaMeshOnProc;

            volMeshOnProc = haveMeshFile
            (
                runTime,
                volMeshMasterInstance/volMeshSubDir,
                "faces"
            );

            if (doFiniteArea)
            {
                areaMeshOnProc = haveMeshFile
                (
                    runTime,
                    areaMeshMasterInstance/areaMeshSubDir,
                    "faceLabels"
                );
            }

            // Prior to loadOrCreateMesh, note which meshes already exist
            // for the current file handler.
            // - where mesh would be written if it didn't exist already.
            fileNameList volMeshDir(Pstream::nProcs());
            {
                volMeshDir[Pstream::myProcNo()] =
                (
                    fileHandler().objectPath
                    (
                        IOobject
                        (
                            "faces",
                            volMeshMasterInstance/volMeshSubDir,
                            runTime
                        ),
                        word::null
                    ).path()
                );

                Pstream::allGatherList(volMeshDir);

                if (optVerbose && Pstream::master())
                {
                    Info<< "Per processor faces dirs:" << nl
                        << '(' << nl;

                    for (const int proci : Pstream::allProcs())
                    {
                        Info<< "    "
                            << runTime.relativePath(volMeshDir[proci]);

                        if (!volMeshOnProc[proci])
                        {
                            Info<< " [missing]";
                        }
                        Info<< nl;
                    }
                    Info<< ')' << nl << endl;
                }
            }

            fileNameList areaMeshDir(Pstream::nProcs());
            if (doFiniteArea)
            {
                areaMeshDir[Pstream::myProcNo()] =
                (
                    fileHandler().objectPath
                    (
                        IOobject
                        (
                            "faceLabels",
                            areaMeshMasterInstance/areaMeshSubDir,
                            runTime
                        ),
                        word::null
                    ).path()
                );

                Pstream::allGatherList(areaMeshDir);

                if (optVerbose && Pstream::master())
                {
                    Info<< "Per processor faceLabels dirs:" << nl
                        << '(' << nl;

                    for (const int proci : Pstream::allProcs())
                    {
                        Info<< "    "
                            << runTime.relativePath(areaMeshDir[proci]);

                        if (!areaMeshOnProc[proci])
                        {
                            Info<< " [missing]";
                        }
                        Info<< nl;
                    }
                    Info<< ')' << nl << endl;
                }
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

            // Volume mesh
            autoPtr<fvMesh> volMeshPtr = fvMeshTools::loadOrCreateMesh
            (
                IOobject
                (
                    regionName,
                    volMeshMasterInstance,
                    runTime,
                    Foam::IOobject::MUST_READ
                ),
                decompose
            );
            fvMesh& mesh = volMeshPtr();


            // Area mesh

            autoPtr<faMesh> areaMeshPtr;

            // Decomposing: must have an undecomposed mesh
            // Redistributing: have any proc mesh
            if
            (
                doFiniteArea
             &&
                (
                    decompose
                  ? areaMeshHaveUndecomposed
                  : areaMeshOnProc.found(true)
                )
            )
            {
                areaMeshPtr = faMeshTools::loadOrCreateMesh
                (
                    IOobject
                    (
                        regionName,
                        areaMeshMasterInstance,
                        runTime,
                        Foam::IOobject::MUST_READ
                    ),
                    mesh,  // <- The referenced polyMesh (from above)
                    decompose
                );

                faMeshTools::forceDemandDriven(*areaMeshPtr);
                faMeshTools::unregisterMesh(*areaMeshPtr);
            }


            if (writeHandler)
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
                           !volMeshOnProc[proci]
                         && volMeshDir[proci] != volMeshDir[myProci]
                        )
                        {
                            Info<< "Deleting mesh dir:"
                                << volMeshDir[proci] << endl;
                            Foam::rmDir(volMeshDir[proci]);
                        }

                        if
                        (
                            !areaMeshOnProc[proci]
                         && areaMeshDir[proci] != areaMeshDir[myProci]
                        )
                        {
                            Info<< "Deleting mesh dir:"
                                << areaMeshDir[proci] << endl;
                            Foam::rmDir(areaMeshDir[proci]);
                        }
                    }

                    // Remove empty directory
                    removeEmptyDirs(mesh.time().path());

                    Pstream::parRun(oldParRun);
                }
            }


            if (Pstream::master() && decompose)
            {
                Info<< "Restoring caseName" << endl;
                runTime.caseName() = proc0CaseName;
                runTime.processorCase(oldProcCase);
            }

            const label nOldCells = mesh.nCells();

            // const label nOldAreaFaces =
            //     (areaMeshPtr ? areaMeshPtr().nFaces() : 0);
            //
            //Pout<< "Loaded mesh : nCells:" << nOldCells
            //    << " nPatches:" << mesh.boundaryMesh().size() << endl;
            //Pout<< "Loaded area mesh : nFaces:" << nOldAreaFaces << endl;

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
                if (!Pstream::master())
                {
                    if (areaMeshPtr && !areaMeshOnProc[Pstream::myProcNo()])
                    {
                        // Remove dummy mesh created by loadOrCreateMesh
                        const bool oldParRun = Pstream::parRun(false);
                        areaMeshPtr->removeFiles();
                        Pstream::parRun(oldParRun);  // Restore parallel state
                    }

                    if (!volMeshOnProc[Pstream::myProcNo()])
                    {
                        // Remove dummy mesh created by loadOrCreateMesh
                        const bool oldParRun = Pstream::parRun(false);
                        mesh.removeFiles();
                        Foam::rmDir(mesh.objectRegistry::objectPath());
                        Pstream::parRun(oldParRun);  // Restore parallel state
                    }
                }
                continue;
            }

            // Area fields first. Read and deregister
            parFaFieldDistributorCache areaFields;
            if (areaMeshPtr)
            {
                areaFields.read
                (
                    baseRunTime,
                    proc0CaseName,
                    decompose,

                    areaMeshOnProc,
                    areaMeshMasterInstance,
                    (*areaMeshPtr)
                );
            }


            // Detect lagrangian fields
            if (Pstream::master() && decompose)
            {
                runTime.caseName() = baseRunTime.caseName();
                runTime.processorCase(false);
            }

            // Read lagrangian fields and store on cloud (objectRegistry)
            PtrList<unmappedPassivePositionParticleCloud> clouds
            (
                readLagrangian
                (
                    mesh,
                    selectedLagrangianFields
                )
            );

            if (Pstream::master() && decompose)
            {
                runTime.caseName() = proc0CaseName;
                runTime.processorCase(oldProcCase);
            }

            // Load fields, do all distribution (mesh and fields)
            // - but not lagrangian fields; these are done later
            autoPtr<mapDistributePolyMesh> distMap = redistributeAndWrite
            (
                std::move(writeHandler),
                baseRunTime,
                proc0CaseName,

                // Controls
                true,           // read fields
                decompose,      // decompose, i.e. read from undecomposed case
                false,          // no reconstruction
                overwrite,

                // Decomposition information
                nDestProcs,
                finalDecomp,

                // For finite volume
                volMeshOnProc,
                volMeshMasterInstance,
                mesh
            );


            // Redistribute any clouds
            redistributeLagrangian
            (
                lagrangianDistributorPtr,
                mesh,
                nOldCells,
                distMap(),
                clouds
            );


            // Redistribute area fields

            mapDistributePolyMesh faDistMap;
            autoPtr<faMesh> areaProcMeshPtr;

            if (areaMeshPtr)
            {
                faDistMap = faMeshDistributor::distribute
                (
                    areaMeshPtr(),
                    distMap(),
                    areaProcMeshPtr
                );

                // Force recreation of everything that might vaguely
                // be used by patches:

                faMeshTools::forceDemandDriven(areaProcMeshPtr());


                if (reconstruct)
                {
                    if (Pstream::master())
                    {
                        Info<< "Setting caseName to " << baseRunTime.caseName()
                            << " to write reconstructed mesh (and fields)."
                            << endl;
                        runTime.caseName() = baseRunTime.caseName();
                        const bool oldProcCase(runTime.processorCase(false));
                        //const bool oldParRun = Pstream::parRun(false);

                        areaProcMeshPtr->write();

                        // Now we've written all. Reset caseName on master
                        Info<< "Restoring caseName" << endl;
                        runTime.caseName() = proc0CaseName;
                        runTime.processorCase(oldProcCase);
                    }
                }
                else
                {
                    autoPtr<fileOperation> defaultHandler;
                    if (writeHandler)
                    {
                        defaultHandler = fileHandler(std::move(writeHandler));
                    }

                    IOmapDistributePolyMeshRef
                    (
                        IOobject
                        (
                            "procAddressing",
                            areaProcMeshPtr->facesInstance(),
                            faMesh::meshSubDir,
                            areaProcMeshPtr->thisDb(),
                            IOobject::NO_READ,
                            IOobject::NO_WRITE,
                            false
                        ),
                        faDistMap
                    ).write();

                    areaProcMeshPtr->write();

                    if (defaultHandler)
                    {
                        writeHandler = fileHandler(std::move(defaultHandler));
                    }

                    if (decompose)
                    {
                        faMeshTools::writeProcAddressing
                        (
                            areaProcMeshPtr(),
                            faDistMap,
                            decompose,
                            std::move(writeHandler)
                        );
                    }
                }

                Info<< "Written redistributed mesh to "
                    << areaProcMeshPtr->facesInstance() << nl << endl;

                faMeshDistributor distributor
                (
                    areaMeshPtr(),      // source
                    areaProcMeshPtr(),  // target
                    faDistMap
                );

                areaFields.redistributeAndWrite(distributor, true);
            }

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
