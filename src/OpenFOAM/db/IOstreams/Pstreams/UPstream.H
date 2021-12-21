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

Class
    Foam::UPstream

Description
    Inter-processor communications stream

SourceFiles
    UPstream.C
    UPstreamCommsStruct.C
    UPstreamTemplates.C
    combineGatherScatter.C
    gatherScatter.C
    gatherScatterList.C

\*---------------------------------------------------------------------------*/

#ifndef UPstream_H
#define UPstream_H

#include "labelList.H"
#include "DynamicList.H"
#include "HashTable.H"
#include "string.H"
#include "Enum.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class UPstream Declaration
\*---------------------------------------------------------------------------*/

class UPstream
{
public:

    //- Int ranges are used for MPI ranks (processes)
    typedef IntRange<int> rangeType;

    //- Types of communications
    enum class commsTypes : char
    {
        blocking,
        scheduled,
        nonBlocking
    };

    //- Names of the communication types
    static const Enum<commsTypes> commsTypeNames;


    // Public Classes

        //- Structure for communicating between processors
        class commsStruct
        {
            // Private Data

                //- procID of above processor
                label above_;

                //- procIDs of processors directly below me
                labelList below_;

                //- procIDs of all processors below (so not just directly below)
                labelList allBelow_;

                //- procIDs of all processors not below.
                //  (inverse set of allBelow_ and minus myProcNo)
                labelList allNotBelow_;


        public:

            // Constructors

                //- Default construct. Above == -1
                commsStruct();

                //- Construct from components
                commsStruct
                (
                    const label above,
                    const labelList& below,
                    const labelList& allBelow,
                    const labelList& allNotBelow
                );

                //- Construct from components; construct allNotBelow_
                commsStruct
                (
                    const label nProcs,
                    const label myProcID,
                    const label above,
                    const labelList& below,
                    const labelList& allBelow
                );


            // Member Functions

                label above() const noexcept
                {
                    return above_;
                }

                const labelList& below() const noexcept
                {
                    return below_;
                }

                const labelList& allBelow() const noexcept
                {
                    return allBelow_;
                }

                const labelList& allNotBelow() const noexcept
                {
                    return allNotBelow_;
                }


            // Member Operators

                bool operator==(const commsStruct&) const;
                bool operator!=(const commsStruct&) const;


             // Ostream Operator

                friend Ostream& operator<<(Ostream&, const commsStruct&);
        };


        //- combineReduce operator for lists. Used for counting.
        struct listEq
        {
            template<class T>
            void operator()(T& x, const T& y) const
            {
                forAll(y, i)
                {
                    if (y[i].size())
                    {
                        x[i] = y[i];
                    }
                }
            }
        };


private:

    // Private Static Data

        //- By default this is not a parallel run
        static bool parRun_;

        //- Have support for threads?
        static bool haveThreads_;

        //- Standard transfer message type
        static int msgType_;

        //- Names of all worlds
        static wordList allWorlds_;

        //- Per processor the name of the world
        static labelList worldIDs_;


        // Communicator specific data

        //- My processor number
        static DynamicList<int> myProcNo_;

        //- List of process IDs
        static DynamicList<List<int>> procIDs_;

        //- Parent communicator
        static DynamicList<label> parentCommunicator_;

        //- Free communicators
        static DynamicList<label> freeComms_;

        //- Linear communication schedule
        static DynamicList<List<commsStruct>> linearCommunication_;

        //- Multi level communication schedule
        static DynamicList<List<commsStruct>> treeCommunication_;


    // Private Member Functions

        //- Set data for parallel running
        static void setParRun(const label nProcs, const bool haveThreads);

        //- Calculate linear communication schedule
        static List<commsStruct> calcLinearComm(const label nProcs);

        //- Calculate tree communication schedule
        static List<commsStruct> calcTreeComm(const label nProcs);

        //- Helper function for tree communication schedule determination
        //  Collects all processorIDs below a processor
        static void collectReceives
        (
            const label procID,
            const List<DynamicList<label>>& receives,
            DynamicList<label>& allReceives
        );

        //- Allocate a communicator with index
        static void allocatePstreamCommunicator
        (
            const label parentIndex,
            const label index
        );

        //- Free a communicator
        static void freePstreamCommunicator
        (
            const label index
        );


protected:

    // Protected Data

        //- Communications type of this stream
        commsTypes commsType_;


public:

    // Declare name of the class and its debug switch
    ClassName("UPstream");


    // Static Data

        //- Should compact transfer be used in which floats replace doubles
        //- reducing the bandwidth requirement at the expense of some loss
        //- in accuracy
        static bool floatTransfer;

        //- Number of processors at which the sum algorithm changes from linear
        //- to tree
        static int nProcsSimpleSum;

        //- Default commsType
        static commsTypes defaultCommsType;

        //- Number of polling cycles in processor updates
        static int nPollProcInterfaces;

        //- Optional maximum message size (bytes)
        static int maxCommsSize;

        //- MPI buffer-size (bytes)
        static const int mpiBufferSize;

        //- Default communicator (all processors)
        static label worldComm;

        //- Debugging: warn for use of any communicator differing from warnComm
        static label warnComm;


    // Constructors

        //- Construct for given communication type
        explicit UPstream(const commsTypes commsType)
        :
            commsType_(commsType)
        {}


    // Member Functions

        //- Allocate a new communicator
        static label allocateCommunicator
        (
            const label parent,
            const labelList& subRanks,
            const bool doPstream = true
        );

        //- Free a previously allocated communicator
        static void freeCommunicator
        (
            const label communicator,
            const bool doPstream = true
        );

        //- Free all communicators
        static void freeCommunicators(const bool doPstream);

        //- Helper class for allocating/freeing communicators
        class communicator
        {
            label comm_;

            //- No copy construct
            communicator(const communicator&) = delete;

            //- No copy assignment
            void operator=(const communicator&) = delete;

        public:

            communicator
            (
                const label parent,
                const labelList& subRanks,
                const bool doPstream
            )
            :
                comm_(allocateCommunicator(parent, subRanks, doPstream))
            {}

            ~communicator()
            {
                freeCommunicator(comm_);
            }

            operator label() const noexcept
            {
                return comm_;
            }
        };

        //- Return physical processor number (i.e. processor number in
        //- worldComm) given communicator and procssor
        static int baseProcNo(const label myComm, const int procID);

        //- Return processor number in communicator (given physical processor
        //- number) (= reverse of baseProcNo)
        static label procNo(const label comm, const int baseProcID);

        //- Return processor number in communicator (given processor number
        //- and communicator)
        static label procNo
        (
            const label myComm,
            const label currentComm,
            const int currentProcID
        );

        //- Add the valid option this type of communications library
        //- adds/requires on the command line
        static void addValidParOptions(HashTable<string>& validParOptions);

        //- Initialisation function called from main
        //  Spawns sub-processes and initialises inter-communication
        static bool init(int& argc, char**& argv, const bool needsThread);

        //- Special purpose initialisation function.
        //  Performs a basic MPI_Init without any other setup.
        //  Only used for applications that need MPI communication when
        //  OpenFOAM is running in a non-parallel mode.
        //  \note Behaves as a no-op if MPI has already been initialized.
        //      Fatal if MPI has already been finalized.
        static bool initNull();


        // Non-blocking comms

            //- Get number of outstanding requests
            static label nRequests();

            //- Truncate number of outstanding requests
            static void resetRequests(const label sz);

            //- Wait until all requests (from start onwards) have finished.
            static void waitRequests(const label start = 0);

            //- Wait until request i has finished.
            static void waitRequest(const label i);

            //- Non-blocking comms: has request i finished?
            static bool finishedRequest(const label i);

            static int allocateTag(const char*);

            static int allocateTag(const word&);

            static void freeTag(const char*, const int tag);

            static void freeTag(const word&, const int tag);


        //- Set as parallel run on/off.
        //  \return the previous value
        static bool parRun(const bool on) noexcept
        {
            bool old(parRun_);
            parRun_ = on;
            return old;
        }

        //- Test if this a parallel run
        //  Modify access is deprecated
        static bool& parRun() noexcept
        {
            return parRun_;
        }

        //- Have support for threads
        static bool haveThreads() noexcept
        {
            return haveThreads_;
        }

        //- Number of processes in parallel run, and 1 for serial run
        static label nProcs(const label communicator = worldComm)
        {
            return procIDs_[communicator].size();
        }

        //- Process index of the master (always 0)
        static constexpr int masterNo() noexcept
        {
            return 0;
        }

        //- Am I the master process
        static bool master(const label communicator = worldComm)
        {
            return myProcNo_[communicator] == masterNo();
        }

        //- Number of this process (starting from masterNo() = 0)
        static int myProcNo(const label communicator = worldComm)
        {
            return myProcNo_[communicator];
        }

        static label parent(const label communicator)
        {
            return parentCommunicator_(communicator);
        }

        //- Process ID of given process index
        static List<int>& procID(label communicator)
        {
            return procIDs_[communicator];
        }


        // Worlds

            //- All worlds
            static const wordList& allWorlds() noexcept
            {
                return allWorlds_;
            }

            //- worldID (index in allWorlds) of all processes
            static const labelList& worldIDs() noexcept
            {
                return worldIDs_;
            }

            //- My worldID
            static label myWorldID()
            {
                return worldIDs_[myProcNo(0)];
            }

            //- My world
            static const word& myWorld()
            {
                return allWorlds()[myWorldID()];
            }


        //- Range of process indices for all processes
        static rangeType allProcs(const label communicator = worldComm)
        {
            // Proc 0 -> nProcs (int value)
            return rangeType(static_cast<int>(nProcs(communicator)));
        }

        //- Range of process indices for sub-processes
        static rangeType subProcs(const label communicator = worldComm)
        {
            // Proc 1 -> nProcs (int value)
            return rangeType(1, static_cast<int>(nProcs(communicator)-1));
        }

        //- Communication schedule for linear all-to-master (proc 0)
        static const List<commsStruct>& linearCommunication
        (
            const label communicator = worldComm
        )
        {
            return linearCommunication_[communicator];
        }

        //- Communication schedule for tree all-to-master (proc 0)
        static const List<commsStruct>& treeCommunication
        (
            const label communicator = worldComm
        )
        {
            return treeCommunication_[communicator];
        }

        //- Message tag of standard messages
        static int& msgType() noexcept
        {
            return msgType_;
        }


        //- Get the communications type of the stream
        commsTypes commsType() const noexcept
        {
            return commsType_;
        }

        //- Set the communications type of the stream
        commsTypes commsType(const commsTypes ct) noexcept
        {
            commsTypes old(commsType_);
            commsType_ = ct;
            return old;
        }


        //- Shutdown (finalize) MPI as required.
        //  Uses MPI_Abort instead of MPI_Finalize if errNo is non-zero
        static void shutdown(int errNo = 0);

        //- Call MPI_Abort with no other checks or cleanup
        static void abort();

        //- Shutdown (finalize) MPI as required and exit program with errNo.
        static void exit(int errNo = 1);

        //- Exchange label with all processors (in the communicator).
        //  sendData[proci] is the label to send to proci.
        //  After return recvData contains the data from the other processors.
        static void allToAll
        (
            const labelUList& sendData,
            labelUList& recvData,
            const label communicator = worldComm
        );

        //- Exchange data with all processors (in the communicator)
        //  sendSizes, sendOffsets give (per processor) the slice of
        //  sendData to send, similarly recvSizes, recvOffsets give the slice
        //  of recvData to receive
        static void allToAll
        (
            const char* sendData,
            const UList<int>& sendSizes,
            const UList<int>& sendOffsets,

            char* recvData,
            const UList<int>& recvSizes,
            const UList<int>& recvOffsets,

            const label communicator = worldComm
        );

        //- Receive data from all processors on the master (low-level)
        static void mpiGather
        (
            const char* sendData,
            int sendSize,

            char* recvData,
            int recvSize,
            const label communicator = worldComm
        );

        //- Send data to all processors from master (low-level)
        static void mpiScatter
        (
            const char* sendData,
            int sendSize,

            char* recvData,
            int recvSize,
            const label communicator = worldComm
        );

        //- Receive data from all processors on the master
        static void gather
        (
            const char* sendData,
            int sendSize,

            char* recvData,
            const UList<int>& recvSizes,
            const UList<int>& recvOffsets,
            const label communicator = worldComm
        );

        //- Send data to all processors from the root of the communicator
        static void scatter
        (
            const char* sendData,
            const UList<int>& sendSizes,
            const UList<int>& sendOffsets,

            char* recvData,
            int recvSize,
            const label communicator = worldComm
        );


    // Gather single, contiguous value(s)

        //- Individual values into list locations.
        //  On master list length == nProcs, otherwise zero length
        template<class T>
        static List<T> listGatherValues
        (
            const T& localValue,
            const label communicator = worldComm
        );

        //- Individual values into list locations.
        //  On master list length == nProcs, otherwise zero length
        template<class T>
        static T listScatterValues
        (
            const UList<T>& allValues,
            const label communicator = worldComm
        );


    // Housekeeping

        //- Process index of first sub-process
        //  \deprecated(2020-09) use subProcs() method instead
        static constexpr int firstSlave() noexcept
        {
            return 1;
        }

        //- Process index of last sub-process
        //  \deprecated(2020-09) use subProcs() method instead
        static int lastSlave(const label communicator = worldComm)
        {
            return nProcs(communicator) - 1;
        }
};


Ostream& operator<<(Ostream&, const UPstream::commsStruct&);

// Template specialisation for access of commsStruct
template<>
UPstream::commsStruct&
UList<UPstream::commsStruct>::operator[](const label);

template<>
const UPstream::commsStruct&
UList<UPstream::commsStruct>::operator[](const label) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "UPstreamTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
