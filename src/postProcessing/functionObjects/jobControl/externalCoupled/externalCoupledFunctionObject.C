/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 OpenCFD Ltd.
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

#include "externalCoupledFunctionObject.H"
#include "addToRunTimeSelectionTable.H"
#include "OSspecific.H"
#include "IFstream.H"
#include "OFstream.H"
#include "volFields.H"
#include "globalIndex.H"
#include "fvMesh.H"
#include "DynamicField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(externalCoupledFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        externalCoupledFunctionObject,
        dictionary
    );
}

Foam::word Foam::externalCoupledFunctionObject::lockName = "OpenFOAM";

Foam::string Foam::externalCoupledFunctionObject::patchKey = "// Patch:";

template<>
const char* Foam::NamedEnum
<
    Foam::externalCoupledFunctionObject::stateEnd,
    2
>::names[] =
{
    "remove",
    "done"
};

const Foam::NamedEnum
<
    Foam::externalCoupledFunctionObject::stateEnd,
    2
> Foam::externalCoupledFunctionObject::stateEndNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::fileName Foam::externalCoupledFunctionObject::baseDir() const
{
    fileName result(commsDir_);
    result.clean();

    return result;
}


Foam::fileName Foam::externalCoupledFunctionObject::groupDir
(
    const fileName& commsDir,
    const word& regionGroupName,
    const wordRe& groupName
)
{
    fileName result
    (
        commsDir
       /regionGroupName
       /string::validate<fileName>(groupName)
    );
    result.clean();

    return result;
}


Foam::fileName Foam::externalCoupledFunctionObject::lockFile() const
{
    return fileName(baseDir()/(lockName + ".lock"));
}


void Foam::externalCoupledFunctionObject::useMaster() const
{
    if (!Pstream::master())
    {
        return;
    }

    const fileName fName(lockFile());
    IFstream is(fName);

    // Only create lock file if it doesn't already exist
    if (!is.good())
    {
        if (log_) Info<< type() << ": creating lock file" << endl;

        OFstream os(fName);
        os  << "status=openfoam\n";
        os.flush();
    }
}


void Foam::externalCoupledFunctionObject::useSlave() const
{
    if (!Pstream::master())
    {
        return;
    }

    if (log_) Info<< type() << ": removing lock file" << endl;

    Foam::rm(lockFile());
}


void Foam::externalCoupledFunctionObject::cleanup() const
{
    if (!Pstream::master())
    {
        return;
    }

    const fileName lck(lockFile());
    switch (stateEnd_)
    {
        case REMOVE:
            {
                if (log_) Info<< type() << ": removing lock file" << endl;
                Foam::rm(lck);
            }
            break;
        case DONE:
            {
                if (log_) Info<< type() << ": lock file status=done" << endl;
                OFstream os(lck);
                os  << "status=done\n";
                os.flush();
            }
            break;
    }
}


void Foam::externalCoupledFunctionObject::removeReadFiles() const
{
    if (!Pstream::master())
    {
        return;
    }

    if (log_) Info<< type() << ": removing all read files" << endl;

    forAll(regionGroupNames_, regionI)
    {
        const word& compName = regionGroupNames_[regionI];

        const labelList& groups = regionToGroups_[compName];
        forAll(groups, i)
        {
            label groupI = groups[i];
            const wordRe& groupName = groupNames_[groupI];

            forAll(groupReadFields_[groupI], fieldI)
            {
                const word& fieldName = groupReadFields_[groupI][fieldI];
                rm
                (
                    groupDir(commsDir_, compName, groupName)
                  / fieldName + ".in"
                );
            }
        }
    }
}


void Foam::externalCoupledFunctionObject::removeWriteFiles() const
{
    if (!Pstream::master())
    {
        return;
    }

    if (log_) Info<< type() << ": removing all write files" << endl;

    forAll(regionGroupNames_, regionI)
    {
        const word& compName = regionGroupNames_[regionI];

        const labelList& groups = regionToGroups_[compName];
        forAll(groups, i)
        {
            label groupI = groups[i];
            const wordRe& groupName = groupNames_[groupI];

            forAll(groupReadFields_[groupI], fieldI)
            {
                const word& fieldName = groupReadFields_[groupI][fieldI];
                rm
                (
                    groupDir(commsDir_, compName, groupName)
                  / fieldName + ".out"
                );
            }
        }
    }
}


void Foam::externalCoupledFunctionObject::waitForSlave() const
{
    const fileName fName(lockFile());
    label totalTime = 0;
    bool found = false;

    if (log_) Info<< type() << ": beginning wait for lock file " << fName << nl;

    while (!found)
    {
        if (Pstream::master())
        {
            if (totalTime > timeOut_)
            {
                FatalErrorInFunction
                    << "Wait time exceeded timeout of " << timeOut_
                    << " s" << abort(FatalError);
            }

            IFstream is(fName);
            if (is.good())
            {
                found = true;

                if (log_)
                {
                    Info<< type() << ": found lock file " << fName << endl;
                }
            }
            else
            {
                sleep(waitInterval_);
                totalTime += waitInterval_;

                if (log_)
                {
                    Info<< type() << ": wait time = " << totalTime << endl;
                }
            }
        }

        // prevent other procs from racing ahead
        reduce(found, orOp<bool>());
    }
}


void Foam::externalCoupledFunctionObject::readColumns
(
    const label nRows,
    const label nColumns,
    autoPtr<IFstream>& masterFilePtr,
    List<scalarField>& data
) const
{
    // Get sizes for all processors
    const globalIndex globalFaces(nRows);

    PstreamBuffers pBufs(Pstream::nonBlocking);
    if (Pstream::master())
    {
        string line;

        // Read data from file and send to destination processor

        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            // Temporary storage
            List<scalarField> values(nColumns);

            // Number of rows to read for processor procI
            label procNRows = globalFaces.localSize(procI);

            forAll(values, columnI)
            {
                values[columnI].setSize(procNRows);
            }

            for (label rowI = 0; rowI < procNRows; rowI++)
            {
                // Get a line
                do
                {
                    if (!masterFilePtr().good())
                    {
                        FatalIOErrorInFunction(masterFilePtr())
                            << "Trying to read data for processor " << procI
                            << " row " << rowI
                            << ". Does your file have as many rows as there are"
                            << " patch faces (" << globalFaces.size()
                            << ") ?" << exit(FatalIOError);
                    }

                    masterFilePtr().getLine(line);
                } while (line.empty() || line[0] == '#');

                IStringStream lineStr(line);

                for (label columnI = 0; columnI < nColumns; columnI++)
                {
                    lineStr >> values[columnI][rowI];
                }
            }

            // Send to procI
            UOPstream str(procI, pBufs);
            str << values;
        }
    }
    pBufs.finishedSends();

    // Read from PstreamBuffers
    UIPstream str(Pstream::masterNo(), pBufs);
    str >> data;
}


void Foam::externalCoupledFunctionObject::readLines
(
    const label nRows,
    autoPtr<IFstream>& masterFilePtr,
    OStringStream& lines
) const
{
    // Get sizes for all processors
    const globalIndex globalFaces(nRows);

    PstreamBuffers pBufs(Pstream::nonBlocking);

    if (Pstream::master())
    {
        string line;

        // Read line from file and send to destination processor

        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            // Number of rows to read for processor procI
            label procNRows = globalFaces.localSize(procI);

            UOPstream toProc(procI, pBufs);

            for (label rowI = 0; rowI < procNRows; rowI++)
            {
                // Get a line
                do
                {
                    if (!masterFilePtr().good())
                    {
                        FatalIOErrorInFunction(masterFilePtr())
                            << "Trying to read data for processor " << procI
                            << " row " << rowI
                            << ". Does your file have as many rows as there are"
                            << " patch faces (" << globalFaces.size()
                            << ") ?" << exit(FatalIOError);
                    }

                    masterFilePtr().getLine(line);
                } while (line.empty() || line[0] == '#');

                // Send line to the destination processor
                toProc << line;
            }
        }
    }


    pBufs.finishedSends();

    // Read lines from PstreamBuffers
    UIPstream str(Pstream::masterNo(), pBufs);
    for (label rowI = 0; rowI < nRows; rowI++)
    {
        string line(str);
        lines << line.c_str() << nl;
    }
}


void Foam::externalCoupledFunctionObject::writeGeometry
(
    const UPtrList<const fvMesh>& meshes,
    const fileName& commsDir,
    const wordRe& groupName
)
{
    wordList regionNames(meshes.size());
    forAll(meshes, i)
    {
        regionNames[i] = meshes[i].dbDir();
    }

    // Make sure meshes are provided in sorted order
    checkOrder(regionNames);

    fileName dir(groupDir(commsDir, compositeName(regionNames), groupName));

    autoPtr<OFstream> osPointsPtr;
    autoPtr<OFstream> osFacesPtr;
    if (Pstream::master())
    {
        mkDir(dir);
        osPointsPtr.reset(new OFstream(dir/"patchPoints"));
        osFacesPtr.reset(new OFstream(dir/"patchFaces"));

        osPointsPtr() << "// Group: " << groupName << endl;
        osFacesPtr()  << "// Group: " << groupName << endl;

        Info<< typeName << ": writing geometry to " << dir << endl;
    }

    // Individual region/patch entries

    DynamicList<face> allFaces;
    DynamicField<point> allPoints;

    labelList pointToGlobal;
    labelList uniquePointIDs;
    forAll(meshes, meshI)
    {
        const fvMesh& mesh = meshes[meshI];

        const labelList patchIDs
        (
            mesh.boundaryMesh().patchSet
            (
                List<wordRe>(1, groupName)
            ).sortedToc()
        );

        forAll(patchIDs, i)
        {
            const polyPatch& p = mesh.boundaryMesh()[patchIDs[i]];

            mesh.globalData().mergePoints
            (
                p.meshPoints(),
                p.meshPointMap(),
                pointToGlobal,
                uniquePointIDs
            );

            label procI = Pstream::myProcNo();

            List<pointField> collectedPoints(Pstream::nProcs());
            collectedPoints[procI] = pointField(mesh.points(), uniquePointIDs);
            Pstream::gatherList(collectedPoints);

            List<faceList> collectedFaces(Pstream::nProcs());
            faceList& patchFaces = collectedFaces[procI];
            patchFaces = p.localFaces();
            forAll(patchFaces, faceI)
            {
                inplaceRenumber(pointToGlobal, patchFaces[faceI]);
            }
            Pstream::gatherList(collectedFaces);

            if (Pstream::master())
            {
                allPoints.clear();
                allFaces.clear();

                for (label procI=0; procI < Pstream::nProcs(); ++procI)
                {
                    allPoints.append(collectedPoints[procI]);
                    allFaces.append(collectedFaces[procI]);
                }

                Info<< typeName << ": mesh " << mesh.name()
                    << ", patch " << p.name()
                    << ": writing " << allPoints.size() << " points to "
                    << osPointsPtr().name() << nl
                    << typeName << ": mesh " << mesh.name()
                    << ", patch " << p.name()
                    << ": writing " << allFaces.size() << " faces to "
                    << osFacesPtr().name() << endl;

                // The entry name (region / patch)
                const string entryHeader =
                    patchKey + ' ' + mesh.name() + ' ' + p.name();

                // Write points
                osPointsPtr()
                    << entryHeader.c_str() << nl << allPoints << nl << endl;

                // Write faces
                osFacesPtr()
                    << entryHeader.c_str() << nl << allFaces << nl << endl;
            }
        }
    }
}


Foam::word Foam::externalCoupledFunctionObject::compositeName
(
    const wordList& regionNames
)
{
    if (regionNames.size() == 0)
    {
        FatalErrorInFunction
            << "Empty regionNames" << abort(FatalError);
        return word::null;
    }
    else if (regionNames.size() == 1)
    {
        if (regionNames[0] == polyMesh::defaultRegion)
        {
            // For compatibility with single region cases suppress single
            // region name
            return word("");
        }
        else
        {
            return regionNames[0];
        }
    }
    else
    {
        // Enforce lexical ordering
        checkOrder(regionNames);

        word composite(regionNames[0]);
        for (label i = 1; i < regionNames.size(); i++)
        {
            composite += "_" + regionNames[i];
        }

        return composite;
    }
}


void Foam::externalCoupledFunctionObject::checkOrder
(
    const wordList& regionNames
)
{
    labelList order;
    sortedOrder(regionNames, order);
    if (order != identity(regionNames.size()))
    {
        FatalErrorInFunction
            << "regionNames " << regionNames << " not in alphabetical order :"
            << order << exit(FatalError);
    }
}


void Foam::externalCoupledFunctionObject::readData()
{
    forAll(regionGroupNames_, regionI)
    {
        const word& compName = regionGroupNames_[regionI];
        const wordList& regionNames = regionGroupRegions_[regionI];

        // Get the meshes for the region-group
        UPtrList<const fvMesh> meshes(regionNames.size());
        forAll(regionNames, j)
        {
            const word& regionName = regionNames[j];
            meshes.set(j, &time_.lookupObject<fvMesh>(regionName));
        }

        const labelList& groups = regionToGroups_[compName];

        forAll(groups, i)
        {
            label groupI = groups[i];
            const wordRe& groupName = groupNames_[groupI];
            const wordList& fieldNames = groupReadFields_[groupI];

            forAll(fieldNames, fieldI)
            {
                const word& fieldName = fieldNames[fieldI];

                bool ok = readData<scalar>
                (
                    meshes,
                    groupName,
                    fieldName
                );
                ok = ok || readData<vector>
                (
                    meshes,
                    groupName,
                    fieldName
                );
                ok = ok || readData<sphericalTensor>
                (
                    meshes,
                    groupName,
                    fieldName
                );
                ok = ok || readData<symmTensor>
                (
                    meshes,
                    groupName,
                    fieldName
                );
                ok = ok || readData<tensor>
                (
                    meshes,
                    groupName,
                    fieldName
                );

                if (!ok)
                {
                    WarningInFunction
                        << "Field " << fieldName << " in regions " << compName
                        << " was not found." << endl;
                }
            }
        }
    }
}


void Foam::externalCoupledFunctionObject::writeData() const
{
    forAll(regionGroupNames_, regionI)
    {
        const word& compName = regionGroupNames_[regionI];
        const wordList& regionNames = regionGroupRegions_[regionI];

        // Get the meshes for the region-group
        UPtrList<const fvMesh> meshes(regionNames.size());
        forAll(regionNames, j)
        {
            const word& regionName = regionNames[j];
            meshes.set(j, &time_.lookupObject<fvMesh>(regionName));
        }

        const labelList& groups = regionToGroups_[compName];

        forAll(groups, i)
        {
            label groupI = groups[i];
            const wordRe& groupName = groupNames_[groupI];
            const wordList& fieldNames = groupWriteFields_[groupI];

            forAll(fieldNames, fieldI)
            {
                const word& fieldName = fieldNames[fieldI];

                bool ok = writeData<scalar>
                (
                    meshes,
                    groupName,
                    fieldName
                );
                ok = ok || writeData<vector>
                (
                    meshes,
                    groupName,
                    fieldName
                );
                ok = ok || writeData<sphericalTensor>
                (
                    meshes,
                    groupName,
                    fieldName
                );
                ok = ok || writeData<symmTensor>
                (
                    meshes,
                    groupName,
                    fieldName
                );
                ok = ok || writeData<tensor>
                (
                    meshes,
                    groupName,
                    fieldName
                );

                if (!ok)
                {
                    WarningInFunction
                        << "Field " << fieldName << " in regions " << compName
                        << " was not found." << endl;
                }
            }
        }
    }
}


void Foam::externalCoupledFunctionObject::initialise()
{
    if (initialised_)
    {
        return;
    }

    // Write the geometry if not already there
    forAll(regionGroupRegions_, i)
    {
        const word& compName = regionGroupNames_[i];
        const wordList& regionNames = regionGroupRegions_[i];

        // Get the meshes for the region-group
        UPtrList<const fvMesh> meshes(regionNames.size());
        forAll(regionNames, j)
        {
            const word& regionName = regionNames[j];
            meshes.set(j, &time_.lookupObject<fvMesh>(regionName));
        }

        const labelList& groups = regionToGroups_[compName];

        forAll(groups, i)
        {
            label groupI = groups[i];
            const wordRe& groupName = groupNames_[groupI];

            bool exists = false;
            if (Pstream::master())
            {
                fileName dir(groupDir(commsDir_, compName, groupName));

                exists =
                    isFile(dir/"patchPoints")
                 || isFile(dir/"patchFaces");
            }

            if (!returnReduce(exists, orOp<bool>()))
            {
                writeGeometry(meshes, commsDir_, groupName);
            }
        }
    }

    if (initByExternal_)
    {
        // Wait for initial data to be made available
        waitForSlave();

        // Read data passed back from external source
        readData();
    }

    initialised_ = true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::externalCoupledFunctionObject::externalCoupledFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObject(name),
    time_(runTime),
    enabled_(true),
    stateEnd_(REMOVE),
    initialised_(false)
{
    read(dict);

    if (Pstream::master())
    {
        mkDir(baseDir());
    }

    if (!initByExternal_)
    {
        useMaster();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::externalCoupledFunctionObject::~externalCoupledFunctionObject()
{
    if (enabled())
    {
        cleanup();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::externalCoupledFunctionObject::on()
{
    enabled_ = true;
}


void Foam::externalCoupledFunctionObject::off()
{
    enabled_ = false;
}


bool Foam::externalCoupledFunctionObject::start()
{
    return true;
}


bool Foam::externalCoupledFunctionObject::execute(const bool forceWrite)
{
    if
    (
        enabled()
     && (!initialised_ || time_.timeIndex() % calcFrequency_ == 0)
    )
    {
        // Initialise the coupling
        initialise();

        // Write data for external source
        writeData();

        // Signal external source to execute (by removing lock file)
        // - Wait for slave to provide data
        useSlave();

        // Wait for response
        waitForSlave();

        // Remove old data files from OpenFOAM
        removeWriteFiles();

        // Read data passed back from external source
        readData();

        // Signal external source to wait (by creating the lock file)
        useMaster();

        return true;
    }
    else
    {
        return false;
    }
}


bool Foam::externalCoupledFunctionObject::end()
{
    if (enabled())
    {
        // Remove old data files
        removeReadFiles();
        removeWriteFiles();
        cleanup();
    }

    return true;
}


bool Foam::externalCoupledFunctionObject::timeSet()
{
    // Do nothing - only valid on execute
    return true;
}


bool Foam::externalCoupledFunctionObject::adjustTimeStep()
{
    return true;
}


bool Foam::externalCoupledFunctionObject::read(const dictionary& dict)
{
    dict.readIfPresent("enabled", enabled_);

    if (!enabled_)
    {
        return true;
    }

    calcFrequency_  = dict.lookupOrDefault("calcFrequency", 1);

    dict.lookup("commsDir") >> commsDir_;
    commsDir_.expand();
    commsDir_.clean();

    waitInterval_   = dict.lookupOrDefault("waitInterval", 1);
    timeOut_        = dict.lookupOrDefault("timeOut", 100*waitInterval_);
    initByExternal_ = readBool(dict.lookup("initByExternal"));
    // initByExternal_ = dict.lookupOrDefault<Switch>("initByExternal", false);
    stateEnd_       =
        stateEndNames_[dict.lookupOrDefault<word>("stateEnd", "remove")];

    log_ = dict.lookupOrDefault("log", false);


    // Get names of all fvMeshes (and derived types)
    wordList allRegionNames(time_.lookupClass<fvMesh>().sortedToc());


    const dictionary& allRegionsDict = dict.subDict("regions");

    forAllConstIter(dictionary, allRegionsDict, iter)
    {
        if (!iter().isDict())
        {
            FatalIOErrorInFunction(allRegionsDict)
                << "Regions must be specified in dictionary format"
                << exit(FatalIOError);
        }

        const wordRe regionGroupName(iter().keyword());
        const dictionary& regionDict = iter().dict();

        labelList regionIDs = findStrings(regionGroupName, allRegionNames);

        const wordList regionNames(allRegionNames, regionIDs);

        regionGroupNames_.append(compositeName(regionNames));
        regionGroupRegions_.append(regionNames);


        forAllConstIter(dictionary, regionDict, regionIter)
        {
            if (!regionIter().isDict())
            {
                FatalIOErrorInFunction(regionDict)
                    << "Regions must be specified in dictionary format"
                    << exit(FatalIOError);
            }
            const wordRe groupName(regionIter().keyword());
            const dictionary& groupDict = regionIter().dict();

            label nGroups = groupNames_.size();
            const wordList readFields(groupDict.lookup("readFields"));
            const wordList writeFields(groupDict.lookup("writeFields"));

            HashTable<labelList>::iterator fnd = regionToGroups_.find
            (
                regionGroupNames_.last()
            );
            if (fnd != regionToGroups_.end())
            {
                fnd().append(nGroups);
            }
            else
            {
                regionToGroups_.insert
                (
                    regionGroupNames_.last(),
                    labelList(1, nGroups)
                );
            }
            groupNames_.append(groupName);
            groupReadFields_.append(readFields);
            groupWriteFields_.append(writeFields);
        }
    }


    // Print a bit
    if (log_)
    {
        Info<< type() << ": Communicating with regions:" << endl;
        forAll(regionGroupNames_, rgI)
        {
            //const wordList& regionNames = regionGroupRegions_[rgI];
            const word& compName = regionGroupNames_[rgI];

            Info<< "Region: " << compName << endl << incrIndent;
            const labelList& groups = regionToGroups_[compName];
            forAll(groups, i)
            {
                label groupI = groups[i];
                const wordRe& groupName = groupNames_[groupI];

                Info<< indent << "patchGroup: " << groupName << "\t"
                    << endl
                    << incrIndent
                    << indent << "Reading fields: "
                    << groupReadFields_[groupI]
                    << endl
                    << indent << "Writing fields: "
                    << groupWriteFields_[groupI]
                    << endl
                    << decrIndent;
            }
            Info<< decrIndent;
        }
        Info<< endl;
    }


    // Note: we should not have to make directories since the geometry
    //       should already be written - but just make sure
    if (Pstream::master())
    {
        forAll(regionGroupNames_, rgI)
        {
            const word& compName = regionGroupNames_[rgI];

            const labelList& groups = regionToGroups_[compName];
            forAll(groups, i)
            {
                label groupI = groups[i];
                const wordRe& groupName = groupNames_[groupI];

                fileName dir(groupDir(commsDir_, compName, groupName));
                if (!isDir(dir))
                {
                    if (log_)
                    {
                        Info<< type() << ": creating communications directory "
                            << dir << endl;
                    }
                    mkDir(dir);
                }
            }
        }
    }

    return true;
}


void Foam::externalCoupledFunctionObject::updateMesh(const mapPolyMesh& mpm)
{}


void Foam::externalCoupledFunctionObject::movePoints(const polyMesh& mesh)
{}


// ************************************************************************* //
