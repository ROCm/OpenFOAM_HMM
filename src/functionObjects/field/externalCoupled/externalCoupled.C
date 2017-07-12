/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2017 OpenCFD Ltd.
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

#include "externalCoupled.H"
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
namespace functionObjects
{
    defineTypeNameAndDebug(externalCoupled, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        externalCoupled,
        dictionary
    );
}
}

const Foam::Enum
<
    Foam::functionObjects::externalCoupled::stateEnd
>
Foam::functionObjects::externalCoupled::stateEndNames_
{
    { stateEnd::REMOVE, "remove" },
    { stateEnd::DONE, "done" }
    // 'IGNORE' is internal use only and thus without a name
};


Foam::word Foam::functionObjects::externalCoupled::lockName = "OpenFOAM";

Foam::string Foam::functionObjects::externalCoupled::patchKey = "// Patch:";


namespace Foam
{
//! \cond fileScope
//- Write list content with size, bracket, content, bracket one-per-line.
//  This makes for consistent for parsing, regardless of the list length.
template <class T>
static void writeList(Ostream& os, const string& header, const UList<T>& L)
{
    // Header string
    os  << header.c_str() << nl;

    // Write size and start delimiter
    os  << L.size() << nl
        << token::BEGIN_LIST << nl;

    // Write contents
    forAll(L, i)
    {
        os << L[i] << nl;
    }

    // Write end delimiter
    os << token::END_LIST << nl << endl;
}
//! \endcond

}
// namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::fileName Foam::functionObjects::externalCoupled::baseDir() const
{
    fileName result(commsDir_);
    result.clean();

    return result;
}


Foam::fileName Foam::functionObjects::externalCoupled::groupDir
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


Foam::fileName Foam::functionObjects::externalCoupled::lockFile() const
{
    return fileName(baseDir()/(lockName + ".lock"));
}


void Foam::functionObjects::externalCoupled::useMaster() const
{
    if (Pstream::master())
    {
        const fileName lck(lockFile());

        // Only create lock file if it doesn't already exist
        if (!Foam::isFile(lck))
        {
            Log << type() << ": creating lock file" << endl;

            OFstream os(lck);
            os << "status=openfoam\n";
            os.flush();
        }
    }
}


void Foam::functionObjects::externalCoupled::useSlave() const
{
    if (Pstream::master())
    {
        Log << type() << ": removing lock file" << endl;

        Foam::rm(lockFile());
    }
}


void Foam::functionObjects::externalCoupled::cleanup() const
{
    if (Pstream::master())
    {
        const fileName lck(lockFile());
        switch (stateEnd_)
        {
            case REMOVE:
            {
                Log << type() << ": removing lock file" << endl;
                Foam::rm(lck);
                break;
            }
            case DONE:
            {
                Log << type() << ": lock file status=done" << endl;
                OFstream os(lck);
                os  << "status=done\n";
                os.flush();
                break;
            }
            case IGNORE:
                break;
        }

        stateEnd_ = IGNORE;   // Avoid re-triggering in destructor
    }
}


void Foam::functionObjects::externalCoupled::removeDataSlave() const
{
    if (!Pstream::master())
    {
        return;
    }

    Log << type() << ": removing data files written by slave" << nl;

    forAll(regionGroupNames_, regioni)
    {
        const word& compName = regionGroupNames_[regioni];

        const labelList& groups = regionToGroups_[compName];
        forAll(groups, i)
        {
            label groupi = groups[i];
            const wordRe& groupName = groupNames_[groupi];

            forAll(groupReadFields_[groupi], fieldi)
            {
                const word& fieldName = groupReadFields_[groupi][fieldi];
                rm
                (
                    groupDir(commsDir_, compName, groupName)
                  / fieldName + ".in"
                );
            }
        }
    }
}


void Foam::functionObjects::externalCoupled::removeDataMaster() const
{
    if (!Pstream::master())
    {
        return;
    }

    Log << type() << ": removing data files written by master" << nl;

    forAll(regionGroupNames_, regioni)
    {
        const word& compName = regionGroupNames_[regioni];

        const labelList& groups = regionToGroups_[compName];
        forAll(groups, i)
        {
            label groupi = groups[i];
            const wordRe& groupName = groupNames_[groupi];

            forAll(groupReadFields_[groupi], fieldi)
            {
                const word& fieldName = groupReadFields_[groupi][fieldi];
                rm
                (
                    groupDir(commsDir_, compName, groupName)
                  / fieldName + ".out"
                );
            }
        }
    }
}


void Foam::functionObjects::externalCoupled::waitForSlave() const
{
    bool waiting = true;
    if (Pstream::master())
    {
        const fileName lck(lockFile());
        unsigned totalTime = 0;

        Log << type() << ": beginning wait for lock file " << lck << nl;

        while ((waiting = !Foam::isFile(lck)) == true)
        {
            sleep(waitInterval_);
            totalTime += waitInterval_;

            if (timeOut_ && totalTime > timeOut_)
            {
                FatalErrorInFunction
                    << "Wait time exceeded timeout of " << timeOut_
                    << " s" << abort(FatalError);
            }

            Log << type() << ": wait time = " << totalTime << endl;
        }

        Log << type() << ": found lock file " << lck << endl;
    }

    // MPI barrier
    Pstream::scatter(waiting);
}


void Foam::functionObjects::externalCoupled::readColumns
(
    const label nRows,
    const label nColumns,
    autoPtr<IFstream>& masterFilePtr,
    List<scalarField>& data
) const
{
    // Get sizes for all processors
    const globalIndex globalFaces(nRows);

    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);
    if (Pstream::master())
    {
        string line;

        // Read data from file and send to destination processor

        for (label proci = 0; proci < Pstream::nProcs(); ++proci)
        {
            // Temporary storage
            List<scalarField> values(nColumns);

            // Number of rows to read for processor proci
            label procNRows = globalFaces.localSize(proci);

            forAll(values, columni)
            {
                values[columni].setSize(procNRows);
            }

            for (label rowi = 0; rowi < procNRows; ++rowi)
            {
                // Get a line
                do
                {
                    if (!masterFilePtr().good())
                    {
                        FatalIOErrorInFunction(masterFilePtr())
                            << "Trying to read data for processor " << proci
                            << " row " << rowi
                            << ". Does your file have as many rows as there are"
                            << " patch faces (" << globalFaces.size()
                            << ") ?" << exit(FatalIOError);
                    }

                    masterFilePtr().getLine(line);
                }
                while (line.empty() || line[0] == '#');

                IStringStream lineStr(line);

                for (label columni = 0; columni < nColumns; ++columni)
                {
                    lineStr >> values[columni][rowi];
                }
            }

            // Send to proci
            UOPstream str(proci, pBufs);
            str << values;
        }
    }
    pBufs.finishedSends();

    // Read from PstreamBuffers
    UIPstream str(Pstream::masterNo(), pBufs);
    str >> data;
}


void Foam::functionObjects::externalCoupled::readLines
(
    const label nRows,
    autoPtr<IFstream>& masterFilePtr,
    OStringStream& lines
) const
{
    // Get sizes for all processors
    const globalIndex globalFaces(nRows);

    PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking);

    if (Pstream::master())
    {
        string line;

        // Read line from file and send to destination processor

        for (label proci = 0; proci < Pstream::nProcs(); ++proci)
        {
            // Number of rows to read for processor proci
            label procNRows = globalFaces.localSize(proci);

            UOPstream toProc(proci, pBufs);

            for (label rowi = 0; rowi < procNRows; ++rowi)
            {
                // Get a line
                do
                {
                    if (!masterFilePtr().good())
                    {
                        FatalIOErrorInFunction(masterFilePtr())
                            << "Trying to read data for processor " << proci
                            << " row " << rowi
                            << ". Does your file have as many rows as there are"
                            << " patch faces (" << globalFaces.size()
                            << ") ?" << exit(FatalIOError);
                    }

                    masterFilePtr().getLine(line);
                }
                while (line.empty() || line[0] == '#');

                // Send line to the destination processor
                toProc << line;
            }
        }
    }


    pBufs.finishedSends();

    // Read lines from PstreamBuffers
    UIPstream str(Pstream::masterNo(), pBufs);
    for (label rowi = 0; rowi < nRows; ++rowi)
    {
        string line(str);
        lines << line.c_str() << nl;
    }
}


void Foam::functionObjects::externalCoupled::writeGeometry
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
    forAll(meshes, meshi)
    {
        const fvMesh& mesh = meshes[meshi];

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

            label proci = Pstream::myProcNo();

            List<pointField> collectedPoints(Pstream::nProcs());
            collectedPoints[proci] = pointField(mesh.points(), uniquePointIDs);
            Pstream::gatherList(collectedPoints);

            List<faceList> collectedFaces(Pstream::nProcs());
            faceList& patchFaces = collectedFaces[proci];
            patchFaces = p.localFaces();
            forAll(patchFaces, facei)
            {
                inplaceRenumber(pointToGlobal, patchFaces[facei]);
            }
            Pstream::gatherList(collectedFaces);

            if (Pstream::master())
            {
                allPoints.clear();
                allFaces.clear();

                for (label proci=0; proci < Pstream::nProcs(); ++proci)
                {
                    allPoints.append(collectedPoints[proci]);
                    allFaces.append(collectedFaces[proci]);
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

                writeList(osPointsPtr(), entryHeader, allPoints);
                writeList(osFacesPtr(),  entryHeader, allFaces);
            }
        }
    }
}


Foam::word Foam::functionObjects::externalCoupled::compositeName
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
            // For compatibility with single region cases
            // - suppress single region name
            return word::null;
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


void Foam::functionObjects::externalCoupled::checkOrder
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


void Foam::functionObjects::externalCoupled::readData()
{
    forAll(regionGroupNames_, regioni)
    {
        const word& compName = regionGroupNames_[regioni];
        const wordList& regionNames = regionGroupRegions_[regioni];

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
            label groupi = groups[i];
            const wordRe& groupName = groupNames_[groupi];
            const wordList& fieldNames = groupReadFields_[groupi];

            forAll(fieldNames, fieldi)
            {
                const word& fieldName = fieldNames[fieldi];

                const bool ok =
                (
                    readData<scalar>
                    (
                        meshes,
                        groupName,
                        fieldName
                    )
                 || readData<vector>
                    (
                        meshes,
                        groupName,
                        fieldName
                    )
                 || readData<sphericalTensor>
                    (
                        meshes,
                        groupName,
                        fieldName
                    )
                 || readData<symmTensor>
                    (
                        meshes,
                        groupName,
                        fieldName
                    )
                 || readData<tensor>
                    (
                        meshes,
                        groupName,
                        fieldName
                    )
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


void Foam::functionObjects::externalCoupled::writeData() const
{
    forAll(regionGroupNames_, regioni)
    {
        const word& compName = regionGroupNames_[regioni];
        const wordList& regionNames = regionGroupRegions_[regioni];

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
            label groupi = groups[i];
            const wordRe& groupName = groupNames_[groupi];
            const wordList& fieldNames = groupWriteFields_[groupi];

            forAll(fieldNames, fieldi)
            {
                const word& fieldName = fieldNames[fieldi];

                const bool ok =
                (
                    writeData<scalar>
                    (
                        meshes,
                        groupName,
                        fieldName
                    )
                 || writeData<vector>
                    (
                        meshes,
                        groupName,
                        fieldName
                    )
                 || writeData<sphericalTensor>
                    (
                        meshes,
                        groupName,
                        fieldName
                    )
                 || writeData<symmTensor>
                    (
                        meshes,
                        groupName,
                        fieldName
                    )
                 || writeData<tensor>
                    (
                        meshes,
                        groupName,
                        fieldName
                    )
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


void Foam::functionObjects::externalCoupled::initialise()
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
            label groupi = groups[i];
            const wordRe& groupName = groupNames_[groupi];

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

    if (slaveFirst_)
    {
        // Wait for initial data to be made available
        waitForSlave();

        // Read data passed back from external source
        readData();
    }

    initialised_ = true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::externalCoupled::externalCoupled
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObject(name),
    time_(runTime),
    stateEnd_(REMOVE),
    initialised_(false)
{
    read(dict);

    if (Pstream::master())
    {
        mkDir(baseDir());
    }

    if (!slaveFirst_)
    {
        useMaster();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::externalCoupled::~externalCoupled()
{
    cleanup();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::externalCoupled::execute()
{
    if (!initialised_ || time_.timeIndex() % calcFrequency_ == 0)
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
        removeDataMaster();

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


bool Foam::functionObjects::externalCoupled::end()
{
    functionObject::end();

    // Remove old data files
    removeDataMaster();
    removeDataSlave();
    cleanup();

    return true;
}


bool Foam::functionObjects::externalCoupled::read(const dictionary& dict)
{
    functionObject::read(dict);

    // NB: Cannot change directory or initialization
    // if things have already been initialized
    if (!initialised_)
    {
        dict.lookup("commsDir") >> commsDir_;
        commsDir_.expand();
        commsDir_.clean();

        slaveFirst_ = readBool(dict.lookup("initByExternal"));
        // slaveFirst_ = dict.lookupOrDefault<bool>("initByExternal", false);
    }

    calcFrequency_ = dict.lookupOrDefault("calcFrequency", 1);

    waitInterval_ = dict.lookupOrDefault("waitInterval", 1u);
    if (!waitInterval_)
    {
        // Enforce non-zero sleep
        waitInterval_ = 1u;
    }

    timeOut_ = dict.lookupOrDefault("timeOut", 100*waitInterval_);
    stateEnd_ =
        stateEndNames_.lookupOrDefault("stateEnd", dict, stateEnd::DONE);


    // Get names of all fvMeshes (and derived types)
    wordList allRegionNames(time_.lookupClass<fvMesh>().sortedToc());

    const dictionary& allRegionsDict = dict.subDict("regions");
    forAllConstIters(allRegionsDict, iter)
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

        forAllConstIters(regionDict, regionIter)
        {
            if (!regionIter().isDict())
            {
                FatalIOErrorInFunction(regionDict)
                    << "Regions must be specified in dictionary format"
                    << exit(FatalIOError);
            }
            const wordRe groupName(regionIter().keyword());
            const dictionary& groupDict = regionIter().dict();

            const label nGroups = groupNames_.size();
            const wordList readFields(groupDict.lookup("readFields"));
            const wordList writeFields(groupDict.lookup("writeFields"));

            auto fnd = regionToGroups_.find(regionGroupNames_.last());
            if (fnd.found())
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


    Info<< type() << ": Communicating with regions:" << endl;
    forAll(regionGroupNames_, rgi)
    {
        //const wordList& regionNames = regionGroupRegions_[rgi];
        const word& compName = regionGroupNames_[rgi];

        Info<< "Region: " << compName << endl << incrIndent;
        const labelList& groups = regionToGroups_[compName];
        forAll(groups, i)
        {
            label groupi = groups[i];
            const wordRe& groupName = groupNames_[groupi];

            Info<< indent << "patchGroup: " << groupName << "\t"
                << endl
                << incrIndent
                << indent << "Reading fields: "
                << groupReadFields_[groupi]
                << endl
                << indent << "Writing fields: "
                << groupWriteFields_[groupi]
                << endl
                << decrIndent;
        }
        Info<< decrIndent;
    }
    Info<< endl;


    // Note: we should not have to make directories since the geometry
    //       should already be written - but just make sure
    if (Pstream::master())
    {
        forAll(regionGroupNames_, rgi)
        {
            const word& compName = regionGroupNames_[rgi];

            const labelList& groups = regionToGroups_[compName];
            forAll(groups, i)
            {
                label groupi = groups[i];
                const wordRe& groupName = groupNames_[groupi];

                fileName dir(groupDir(commsDir_, compName, groupName));
                if (!isDir(dir))
                {
                    Log << type() << ": creating communications directory "
                        << dir << endl;
                    mkDir(dir);
                }
            }
        }
    }

    return true;
}


bool Foam::functionObjects::externalCoupled::write()
{
    return true;
}


// ************************************************************************* //
