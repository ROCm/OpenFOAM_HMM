/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2018 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

#include "collatedFileOperation.H"
#include "addToRunTimeSelectionTable.H"
#include "Pstream.H"
#include "Time.H"
#include "threadedCollatedOFstream.H"
#include "decomposedBlockData.H"
#include "registerSwitch.H"
#include "masterOFstream.H"
#include "OFstream.H"
#include "foamVersion.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
namespace fileOperations
{
    defineTypeNameAndDebug(collatedFileOperation, 0);
    addToRunTimeSelectionTable
    (
        fileOperation,
        collatedFileOperation,
        word
    );

    float collatedFileOperation::maxThreadFileBufferSize
    (
        debug::floatOptimisationSwitch("maxThreadFileBufferSize", 1e9)
    );
    registerOptSwitch
    (
        "maxThreadFileBufferSize",
        float,
        collatedFileOperation::maxThreadFileBufferSize
    );

    // Mark as needing threaded mpi
    addNamedToRunTimeSelectionTable
    (
        fileOperationInitialise,
        collatedFileOperationInitialise,
        word,
        collated
    );
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fileOperations::collatedFileOperation::printBanner
(
    const bool printRanks
) const
{
    DetailInfo
        << "I/O    : " << this->type();

    if (maxThreadFileBufferSize == 0)
    {
        DetailInfo
            << " [unthreaded] (maxThreadFileBufferSize = 0)." << nl
            << "         Writing may be slow for large file sizes."
            << endl;
    }
    else
    {
        DetailInfo
            << " [threaded] (maxThreadFileBufferSize = "
            << maxThreadFileBufferSize << ")." << nl
            << "         Requires buffer large enough to collect all data"
               " or thread support" << nl
            << "         enabled in MPI. If MPI thread support cannot be"
               " enabled, deactivate" << nl
            << "         threading by setting maxThreadFileBufferSize"
               " to 0 in" << nl
            << "         OpenFOAM etc/controlDict" << endl;
    }

    if (printRanks)
    {
        // Information about the ranks
        stringList hosts(Pstream::nProcs());
        if (Pstream::master(comm_))
        {
            // Don't usually need the pid
            // hosts[Pstream::myProcNo()] = hostName()+"."+name(pid());
            hosts[Pstream::myProcNo()] = hostName();
        }
        Pstream::gatherList(hosts);

        DynamicList<label> offsetMaster(Pstream::nProcs());

        forAll(hosts, ranki)
        {
            if (!hosts[ranki].empty())
            {
                offsetMaster.append(ranki);
            }
        }

        if (offsetMaster.size() > 1)
        {
            DetailInfo
                << "IO nodes:" << nl << '(' << nl;

            offsetMaster.append(Pstream::nProcs());

            for (label group = 1; group < offsetMaster.size(); ++group)
            {
                const label beg = offsetMaster[group-1];
                const label end = offsetMaster[group];

                DetailInfo
                    << "    (" << hosts[beg].c_str() << ' '
                    << (end-beg) << ')' << nl;
            }
            DetailInfo
                << ')' << nl;
        }
    }

    // if (IOobject::fileModificationChecking == IOobject::timeStampMaster)
    // {
    //     WarningInFunction
    //         << "Resetting fileModificationChecking to timeStamp" << endl;
    // }
    // else if (IOobject::fileModificationChecking == IOobject::inotifyMaster)
    // {
    //     WarningInFunction
    //         << "Resetting fileModificationChecking to inotify" << endl;
    // }
}


bool Foam::fileOperations::collatedFileOperation::isMasterRank
(
    const label proci
)
const
{
    if (Pstream::parRun())
    {
        return Pstream::master(comm_);
    }
    else if (ioRanks_.size())
    {
        // Found myself in IO rank
        return ioRanks_.found(proci);
    }
    else
    {
        // Assume all in single communicator
        return proci == 0;
    }
}


bool Foam::fileOperations::collatedFileOperation::appendObject
(
    const regIOobject& io,
    const fileName& pathName,
    IOstreamOption streamOpt
) const
{
    // Append to processorsNN/ file

    const label proci = detectProcessorPath(io.objectPath());

    if (debug)
    {
        Pout<< "collatedFileOperation::writeObject :"
            << " For local object : " << io.name()
            << " appending processor " << proci
            << " data to " << pathName << endl;
    }
    if (proci == -1)
    {
        FatalErrorInFunction
            << "Invalid processor path: " << pathName
            << exit(FatalError);
    }

    const bool isMaster = isMasterRank(proci);

    // Update meta-data for current state
    if (isMaster)
    {
        const_cast<regIOobject&>(io).updateMetaData();
    }

    // Note: cannot do append + compression. This is a limitation
    // of ogzstream (or rather most compressed formats)

    OFstream os
    (
        pathName,
        IOstreamOption(IOstream::BINARY, streamOpt.version()),  // UNCOMPRESSED
        !isMaster  // append slaves
    );

    if (!os.good())
    {
        FatalIOErrorInFunction(os)
            << "Cannot open for appending"
            << exit(FatalIOError);
    }

    if (isMaster)
    {
        decomposedBlockData::writeHeader(os, streamOpt, io);
    }

    std::streamoff blockOffset = decomposedBlockData::writeBlockEntry
    (
        os,
        streamOpt,
        io,
        proci,
        // With FoamFile header on master?
        isMaster
    );

    return (blockOffset >= 0) && os.good();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileOperations::collatedFileOperation::collatedFileOperation
(
    bool verbose
)
:
    masterUncollatedFileOperation
    (
        (
            ioRanks().size()
          ? UPstream::allocateCommunicator
            (
                UPstream::worldComm,
                subRanks(Pstream::nProcs())
            )
          : UPstream::worldComm
        ),
        false
    ),
    myComm_(comm_),
    writer_(mag(maxThreadFileBufferSize), comm_),
    nProcs_(Pstream::nProcs()),
    ioRanks_(ioRanks())
{
    if (verbose && Foam::infoDetailLevel > 0)
    {
        this->printBanner(ioRanks_.size());
    }
}


Foam::fileOperations::collatedFileOperation::collatedFileOperation
(
    const label comm,
    const labelList& ioRanks,
    const word& typeName,
    bool verbose
)
:
    masterUncollatedFileOperation(comm, false),
    myComm_(-1),
    writer_(mag(maxThreadFileBufferSize), comm),
    nProcs_(Pstream::nProcs()),
    ioRanks_(ioRanks)
{
    if (verbose && Foam::infoDetailLevel > 0)
    {
        this->printBanner(ioRanks_.size());
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fileOperations::collatedFileOperation::~collatedFileOperation()
{
    // Wait for any outstanding file operations
    flush();

    if (myComm_ != -1 && myComm_ != UPstream::worldComm)
    {
        UPstream::freeCommunicator(myComm_);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::fileOperations::collatedFileOperation::objectPath
(
    const IOobject& io,
    const word& typeName
) const
{
    // Replacement for objectPath
    if (io.time().processorCase())
    {
        return masterUncollatedFileOperation::localObjectPath
        (
            io,
            fileOperation::PROCOBJECT,
            "dummy",        // not used for processorsobject
            io.instance()
        );
    }
    else
    {
        return masterUncollatedFileOperation::localObjectPath
        (
            io,
            fileOperation::OBJECT,
            word::null,
            io.instance()
        );
    }
}


bool Foam::fileOperations::collatedFileOperation::writeObject
(
    const regIOobject& io,
    IOstreamOption streamOpt,
    const bool valid
) const
{
    const Time& tm = io.time();
    const fileName& inst = io.instance();

    // Update meta-data for current state
    const_cast<regIOobject&>(io).updateMetaData();

    if (inst.isAbsolute() || !tm.processorCase())
    {
        mkDir(io.path());
        fileName pathName(io.objectPath());

        if (debug)
        {
            Pout<< "collatedFileOperation::writeObject :"
                << " For object : " << io.name()
                << " falling back to master-only output to " << io.path()
                << endl;
        }

        masterOFstream os
        (
            pathName,
            streamOpt,
            false,  // append=false
            valid
        );

        // If any of these fail, return
        // (leave error handling to Ostream class)

        const bool ok =
        (
            os.good()
         && io.writeHeader(os)
         && io.writeData(os)
        );

        if (ok)
        {
            IOobject::writeEndDivider(os);
        }

        return ok;
    }
    else
    {
        // Construct the equivalent processors/ directory
        fileName path(processorsPath(io, inst, processorsDir(io)));

        mkDir(path);
        fileName pathName(path/io.name());

        if (io.global())
        {
            if (debug)
            {
                Pout<< "collatedFileOperation::writeObject :"
                    << " For global object : " << io.name()
                    << " falling back to master-only output to " << pathName
                    << endl;
            }

            masterOFstream os
            (
                pathName,
                streamOpt,
                false,  // append=false
                valid
            );

            // If any of these fail, return
            // (leave error handling to Ostream class)

            const bool ok =
            (
                os.good()
             && io.writeHeader(os)
             && io.writeData(os)
            );

            if (ok)
            {
                IOobject::writeEndDivider(os);
            }

            return ok;
        }
        else if (!Pstream::parRun())
        {
            // Special path for e.g. decomposePar. Append to
            // processorsDDD/ file
            if (debug)
            {
                Pout<< "collatedFileOperation::writeObject :"
                    << " For object : " << io.name()
                    << " appending to " << pathName << endl;
            }

            return appendObject(io, pathName, streamOpt);
        }
        else
        {
            // Re-check static maxThreadFileBufferSize variable to see
            // if needs to use threading
            const bool useThread = (maxThreadFileBufferSize != 0);

            if (debug)
            {
                Pout<< "collatedFileOperation::writeObject :"
                    << " For object : " << io.name()
                    << " starting collating output to " << pathName
                    << " useThread:" << useThread << endl;
            }

            if (!useThread)
            {
                writer_.waitAll();
            }

            threadedCollatedOFstream os
            (
                writer_,
                pathName,
                streamOpt,
                useThread
            );

            bool ok = os.good();

            if (Pstream::master(comm_))
            {
                // Suppress comment banner
                const bool old = IOobject::bannerEnabled(false);

                ok = ok && io.writeHeader(os);

                IOobject::bannerEnabled(old);

                // Additional header content
                dictionary dict;
                decomposedBlockData::writeExtraHeaderContent
                (
                    dict,
                    streamOpt,
                    io
                );
                os.setHeaderEntries(dict);
            }

            ok = ok && io.writeData(os);
            // No end divider for collated output

            return ok;
        }
    }
}

void Foam::fileOperations::collatedFileOperation::flush() const
{
    if (debug)
    {
        Pout<< "collatedFileOperation::flush : clearing and waiting for thread"
            << endl;
    }
    masterUncollatedFileOperation::flush();
    // Wait for thread to finish (note: also removes thread)
    writer_.waitAll();
}


Foam::word Foam::fileOperations::collatedFileOperation::processorsDir
(
    const fileName& fName
) const
{
    if (Pstream::parRun())
    {
        const List<int>& procs(UPstream::procID(comm_));

        word procDir(processorsBaseDir+Foam::name(Pstream::nProcs()));

        if (procs.size() != Pstream::nProcs())
        {
            procDir +=
              + "_"
              + Foam::name(procs.first())
              + "-"
              + Foam::name(procs.last());
        }
        return procDir;
    }
    else
    {
        word procDir(processorsBaseDir+Foam::name(nProcs_));

        if (ioRanks_.size())
        {
            // Detect current processor number
            label proci = detectProcessorPath(fName);

            if (proci != -1)
            {
                // Find lowest io rank
                label minProc = 0;
                label maxProc = nProcs_-1;
                for (const label ranki : ioRanks_)
                {
                    if (ranki >= nProcs_)
                    {
                        break;
                    }
                    else if (ranki <= proci)
                    {
                        minProc = ranki;
                    }
                    else
                    {
                        maxProc = ranki-1;
                        break;
                    }
                }
                procDir +=
                  + "_"
                  + Foam::name(minProc)
                  + "-"
                  + Foam::name(maxProc);
            }
        }

        return procDir;
    }
}


Foam::word Foam::fileOperations::collatedFileOperation::processorsDir
(
    const IOobject& io
) const
{
    return processorsDir(io.objectPath());
}


void Foam::fileOperations::collatedFileOperation::setNProcs(const label nProcs)
{
    nProcs_ = nProcs;

    if (debug)
    {
        Pout<< "collatedFileOperation::setNProcs :"
            << " Setting number of processors to " << nProcs_ << endl;
    }
}


// ************************************************************************* //
