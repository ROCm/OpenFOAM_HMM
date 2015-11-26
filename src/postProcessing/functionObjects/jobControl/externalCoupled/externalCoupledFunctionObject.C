/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenCFD Ltd.
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

Foam::string Foam::externalCoupledFunctionObject::patchKey = "# Patch: ";


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
    const word& regionName,
    const wordRe& groupName
)
{
    fileName result(commsDir/regionName/string::validate<fileName>(groupName));
    result.clean();

    return result;
}


Foam::fileName Foam::externalCoupledFunctionObject::lockFile() const
{
    return fileName(baseDir()/(lockName + ".lock"));
}


void Foam::externalCoupledFunctionObject::createLockFile() const
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
        os  << "lock file";
        os.flush();
    }
}


void Foam::externalCoupledFunctionObject::removeLockFile() const
{
    if (!Pstream::master())
    {
        return;
    }

    if (log_) Info<< type() << ": removing lock file" << endl;

    rm(lockFile());
}


void Foam::externalCoupledFunctionObject::removeReadFiles() const
{
    if (!Pstream::master())
    {
        return;
    }

    if (log_) Info<< type() << ": removing all read files" << endl;

    forAll(regionNames_, regionI)
    {
        const word& regionName = regionNames_[regionI];
        const fvMesh& mesh = time_.lookupObject<fvMesh>(regionName);
        const labelList& groups = regionToGroups_[regionName];
        forAll(groups, i)
        {
            label groupI = groups[i];
            const wordRe& groupName = groupNames_[groupI];

            forAll(groupReadFields_[groupI], fieldI)
            {
                const word& fieldName = groupReadFields_[groupI][fieldI];
                rm
                (
                    groupDir(commsDir_, mesh.dbDir(), groupName)
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

    forAll(regionNames_, regionI)
    {
        const word& regionName = regionNames_[regionI];
        const fvMesh& mesh = time_.lookupObject<fvMesh>(regionName);
        const labelList& groups = regionToGroups_[regionName];
        forAll(groups, i)
        {
            label groupI = groups[i];
            const wordRe& groupName = groupNames_[groupI];

            forAll(groupWriteFields_[groupI], fieldI)
            {
                const word& fieldName = groupWriteFields_[groupI][fieldI];
                rm
                (
                    groupDir(commsDir_, mesh.dbDir(), groupName)
                  / fieldName + ".out"
                );
            }
        }
    }
}


void Foam::externalCoupledFunctionObject::wait() const
{
    const fileName fName(lockFile());
    label found = 0;
    label totalTime = 0;

    if (log_) Info<< type() << ": beginning wait for lock file " << fName << nl;

    while (found == 0)
    {
        if (Pstream::master())
        {
            if (totalTime > timeOut_)
            {
                FatalErrorIn
                (
                    "void "
                    "Foam::externalCoupledFunctionObject::wait() "
                    "const"
                )
                    << "Wait time exceeded time out time of " << timeOut_
                    << " s" << abort(FatalError);
            }

            IFstream is(fName);

            if (is.good())
            {
                found++;

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
        reduce(found, sumOp<label>());
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
                        FatalIOErrorIn
                        (
                            "externalCoupledFunctionObject::readColumns()",
                            masterFilePtr()
                        )   << "Trying to read data for processor " << procI
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
                        FatalIOErrorIn
                        (
                            "externalCoupledFunctionObject::readColumns()",
                            masterFilePtr()
                        )   << "Trying to read data for processor " << procI
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
    const fvMesh& mesh,
    const fileName& commsDir,
    const wordRe& groupName
)
{
    fileName dir(groupDir(commsDir, mesh.dbDir(), groupName));

    //if (log_)
    {
        Info<< typeName << ": writing geometry to " << dir << endl;
    }

    autoPtr<OFstream> osPointsPtr;
    autoPtr<OFstream> osFacesPtr;
    if (Pstream::master())
    {
        mkDir(dir);
        osPointsPtr.reset(new OFstream(dir/"patchPoints"));
        osFacesPtr.reset(new OFstream(dir/"patchFaces"));
    }

    const labelList patchIDs
    (
        mesh.boundaryMesh().patchSet
        (
            List<wordRe>(1, groupName)
        ).sortedToc()
    );

    forAll(patchIDs, i)
    {
        label patchI = patchIDs[i];

        const polyPatch& p = mesh.boundaryMesh()[patchI];

        labelList pointToGlobal;
        labelList uniquePointIDs;
        mesh.globalData().mergePoints
        (
            p.meshPoints(),
            p.meshPointMap(),
            pointToGlobal,
            uniquePointIDs
        );

        label procI = Pstream::myProcNo();

        List<pointField> allPoints(Pstream::nProcs());
        allPoints[procI] = pointField(mesh.points(), uniquePointIDs);
        Pstream::gatherList(allPoints);

        List<faceList> allFaces(Pstream::nProcs());
        faceList& patchFaces = allFaces[procI];
        patchFaces = p.localFaces();
        forAll(patchFaces, faceI)
        {
            inplaceRenumber(pointToGlobal, patchFaces[faceI]);
        }
        Pstream::gatherList(allFaces);

        if (Pstream::master())
        {
            pointField pts
            (
                ListListOps::combine<pointField>
                (
                    allPoints,
                    accessOp<pointField>()
                )
            );

            //if (log_)
            {
                Info<< typeName << ": for patch " << p.name()
                    << " writing " << pts.size() << " points to "
                    << osPointsPtr().name() << endl;
            }

            // Write points
            osPointsPtr() << patchKey.c_str() << p.name() << pts << endl;

            faceList fcs
            (
                ListListOps::combine<faceList>(allFaces, accessOp<faceList>())
            );

            //if (log_)
            {
                Info<< typeName << ": for patch " << p.name()
                    << " writing " << fcs.size() << " faces to "
                    << osFacesPtr().name() << endl;
            }

            // Write faces
            osFacesPtr() << patchKey.c_str() << p.name() << fcs << endl;
        }
    }
}


void Foam::externalCoupledFunctionObject::readData()
{
    forAll(regionNames_, regionI)
    {
        const word& regionName = regionNames_[regionI];
        const labelList& groups = regionToGroups_[regionName];

        const fvMesh& mesh = time_.lookupObject<fvMesh>(regionName);

        forAll(groups, i)
        {
            label groupI = groups[i];
            const wordRe& groupName = groupNames_[groupI];
            const labelList& patchIDs = groupPatchIDs_[groupI];
            const wordList& fieldNames = groupReadFields_[groupI];

            forAll(fieldNames, fieldI)
            {
                const word& fieldName = fieldNames[fieldI];

                bool ok = readData<scalar>
                (
                    mesh,
                    groupName,
                    patchIDs,
                    fieldName
                );
                ok = ok || readData<vector>
                (
                    mesh,
                    groupName,
                    patchIDs,
                    fieldName
                );
                ok = ok || readData<sphericalTensor>
                (
                    mesh,
                    groupName,
                    patchIDs,
                    fieldName
                );
                ok = ok || readData<symmTensor>
                (
                    mesh,
                    groupName,
                    patchIDs,
                    fieldName
                );
                ok = ok || readData<tensor>
                (
                    mesh,
                    groupName,
                    patchIDs,
                    fieldName
                );

                if (!ok)
                {
                    WarningIn
                    (
                        "void Foam::externalCoupledFunctionObject::readData()"
                    )
                        << "Field " << fieldName << " in region " << mesh.name()
                        << " was not found." << endl;
                }
            }
        }
    }
}


void Foam::externalCoupledFunctionObject::writeData() const
{
    forAll(regionNames_, regionI)
    {
        const word& regionName = regionNames_[regionI];
        const labelList& groups = regionToGroups_[regionName];

        const fvMesh& mesh = time_.lookupObject<fvMesh>(regionName);

        forAll(groups, i)
        {
            label groupI = groups[i];
            const wordRe& groupName = groupNames_[groupI];
            const labelList& patchIDs = groupPatchIDs_[groupI];
            const wordList& fieldNames = groupWriteFields_[groupI];

            forAll(fieldNames, fieldI)
            {
                const word& fieldName = fieldNames[fieldI];
                bool ok = writeData<scalar>
                (
                    mesh,
                    groupName,
                    patchIDs,
                    fieldName
                );
                ok = ok || writeData<vector>
                (
                    mesh,
                    groupName,
                    patchIDs,
                    fieldName
                );
                ok = ok || writeData<sphericalTensor>
                (
                    mesh,
                    groupName,
                    patchIDs,
                    fieldName
                );
                ok = ok || writeData<symmTensor>
                (
                    mesh,
                    groupName,
                    patchIDs,
                    fieldName
                );
                ok = ok || writeData<tensor>
                (
                    mesh,
                    groupName,
                    patchIDs,
                    fieldName
                );

                if (!ok)
                {
                    WarningIn
                    (
                        "void Foam::externalCoupledFunctionObject::writeData()"
                    )
                        << "Field " << fieldName << " in region " << mesh.name()
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
    forAll(regionNames_, regionI)
    {
        const word& regionName = regionNames_[regionI];
        const labelList& groups = regionToGroups_[regionName];

        const fvMesh& mesh = time_.lookupObject<fvMesh>(regionName);

        forAll(groups, i)
        {
            label groupI = groups[i];
            const wordRe& groupName = groupNames_[groupI];

            bool exists = false;
            if (Pstream::master())
            {
                fileName dir(groupDir(commsDir_, mesh.dbDir(), groupName));

                exists = isFile(dir/"patchPoints") || isFile(dir/"patchFaces");
            }

            if (!returnReduce(exists, orOp<bool>()))
            {
                writeGeometry(mesh, commsDir_, groupName);
            }
        }
    }

    if (initByExternal_)
    {
        // Wait for initial data to be made available
        wait();

        // Eead data passed back from external source
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
    initialised_(false)
{
    read(dict);

    if (Pstream::master())
    {
        mkDir(baseDir());
    }

    if (!initByExternal_)
    {
        createLockFile();
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::externalCoupledFunctionObject::~externalCoupledFunctionObject()
{}


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

        // remove lock file, signalling external source to execute
        removeLockFile();

        // Wait for response
        wait();

        // Remove old data files from OpenFOAM
        removeWriteFiles();

        // Read data passed back from external source
        readData();

        // create lock file for external source
        createLockFile();

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
        removeLockFile();
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

    dict.lookup("commsDir") >> commsDir_;
    commsDir_.expand();

    waitInterval_ = dict.lookupOrDefault("waitInterval", 1);
    timeOut_ = dict.lookupOrDefault("timeOut", 100*waitInterval_);
    calcFrequency_ = dict.lookupOrDefault("calcFrequency", 1);
    initByExternal_ = readBool(dict.lookup("initByExternal"));
    log_ = dict.lookupOrDefault("log", false);

    const dictionary& allRegionsDict = dict.subDict("regions");

    forAllConstIter(dictionary, allRegionsDict, iter)
    {
        if (!iter().isDict())
        {
            FatalIOErrorIn
            (
                "void Foam::externalCoupledFunctionObject::read"
                "(const dictionary&)",
                allRegionsDict
            )
                << "Regions must be specified in dictionary format"
                << exit(FatalIOError);
        }

        const word& regionName = iter().keyword();
        const dictionary& regionDict = iter().dict();
        regionNames_.append(regionName);

        forAllConstIter(dictionary, regionDict, regionIter)
        {
            if (!regionIter().isDict())
            {
                FatalIOErrorIn
                (
                    "void Foam::externalCoupledFunctionObject::read"
                    "(const dictionary&)",
                    regionDict
                )
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
                regionName
            );
            if (fnd != regionToGroups_.end())
            {
                fnd().append(nGroups);
            }
            else
            {
                regionToGroups_.insert(regionName, labelList(1, nGroups));
            }
            groupNames_.append(groupName);
            groupReadFields_.append(readFields);
            groupWriteFields_.append(writeFields);

            // Pre-calculate the patchIDs
            const fvMesh& mesh = time_.lookupObject<fvMesh>(regionName);
            groupPatchIDs_.append
            (
                mesh.boundaryMesh().patchSet
                (
                    List<wordRe>(1, groupName)
                ).sortedToc()
            );
        }
    }


    // Print a bit
    if (log_)
    {
        Info<< type() << ": Communicating with regions:" << endl;
        forAll(regionNames_, regionI)
        {
            const word& regionName = regionNames_[regionI];
            const fvMesh& mesh = time_.lookupObject<fvMesh>(regionName);

            Info<< "Region: " << mesh.name() << endl << incrIndent;
            const labelList& groups = regionToGroups_[regionName];
            forAll(groups, i)
            {
                label groupI = groups[i];
                const wordRe& groupName = groupNames_[groupI];
                const labelList& patchIDs = groupPatchIDs_[groupI];

                Info<< indent << "Group: " << groupName << "\t"
                    << " patches: " << patchIDs << endl
                    << incrIndent
                    << indent << "Reading fields: " << groupReadFields_[groupI]
                    << endl
                    << indent << "Writing fields: " << groupWriteFields_[groupI]
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
        forAll(regionNames_, regionI)
        {
            const word& regionName = regionNames_[regionI];
            const fvMesh& mesh = time_.lookupObject<fvMesh>(regionName);
            const labelList& groups = regionToGroups_[regionName];
            forAll(groups, i)
            {
                label groupI = groups[i];
                const wordRe& groupName = groupNames_[groupI];

                fileName dir(groupDir(commsDir_, mesh.dbDir(), groupName));
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
