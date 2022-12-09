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

Class
    Foam::UPstream

Description
    Inter-processor communications stream

SourceFiles
    UPstream.C
    UPstreamCommsStruct.C
    UPstreamTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_UPstream_H
#define Foam_UPstream_H

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
        blocking,       //!< "blocking"    : (MPI_Bsend, MPI_Recv)
        scheduled,      //!< "scheduled"   : (MPI_Send, MPI_Recv)
        nonBlocking     //!< "nonBlocking" : (MPI_Isend, MPI_Irecv)
    };

    //- Enumerated names for the communication types
    static const Enum<commsTypes> commsTypeNames;


    // Public Classes

        //- Structure for communicating between processors
        class commsStruct
        {
            // Private Data

                //- The procID of the processor directly above
                label above_;

                //- The procIDs of all processors directly below
                labelList below_;

                //- procIDs of all processors below (so not just directly below)
                labelList allBelow_;

                //- procIDs of all processors not below.
                //  Inverse set of allBelow_ without myProcNo.
                labelList allNotBelow_;


        public:

            // Constructors

                //- Default construct. Above == -1
                commsStruct() noexcept;

                //- Construct from components
                commsStruct
                (
                    const label above,
                    const labelUList& below,
                    const labelUList& allBelow,
                    const labelUList& allNotBelow
                );

                //- Construct from components; construct allNotBelow_
                commsStruct
                (
                    const label nProcs,
                    const label myProcID,
                    const label above,
                    const labelUList& below,
                    const labelUList& allBelow
                );


            // Member Functions

                //- Reset to default constructed state
                void clear();

                //- The procID of the processor directly above
                label above() const noexcept
                {
                    return above_;
                }

                //- The procIDs of all processors directly below
                const labelList& below() const noexcept
                {
                    return below_;
                }

                //- The procIDs of all processors below
                //- (so not just directly below)
                const labelList& allBelow() const noexcept
                {
                    return allBelow_;
                }

                //- The procIDs of all processors that are above.
                //- The inverse set of allBelow without myProcNo.
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

    // Private Data

        //- Communications type of this stream
        commsTypes commsType_;


    // Private Static Data

        //- By default this is not a parallel run
        static bool parRun_;

        //- Have support for threads?
        static bool haveThreads_;

        //- Standard transfer message type
        static int msgType_;

        //- Names of all worlds
        static wordList allWorlds_;

        //- Per processor the world index (into allWorlds_)
        static labelList worldIDs_;


    // Communicator specific data

        //- My processor number
        static DynamicList<int> myProcNo_;

        //- List of process IDs
        static DynamicList<List<int>> procIDs_;

        //- Parent communicator
        static DynamicList<label> parentComm_;

        //- Free communicators
        static DynamicList<label> freeComms_;

        //- Linear communication schedule
        static DynamicList<List<commsStruct>> linearCommunication_;

        //- Multi level communication schedule
        static DynamicList<List<commsStruct>> treeCommunication_;


    // Private Member Functions

        //- Set data for parallel running
        static void setParRun(const label nProcs, const bool haveThreads);

        //- Allocate a communicator with index
        static void allocatePstreamCommunicator
        (
            const label parentIndex,
            const label index
        );

        //- Free MPI components of communicator.
        //  Does not touch the first two communicators (SELF, WORLD)
        static void freePstreamCommunicator(const label index);


public:

    //- Declare name of the class and its debug switch
    ClassName("UPstream");


    // Static Data

        //- Should compact transfer be used in which floats replace doubles
        //- reducing the bandwidth requirement at the expense of some loss
        //- in accuracy
        static bool floatTransfer;

        //- Number of processors to change from linear to tree communication
        static int nProcsSimpleSum;

        //- Default commsType
        static commsTypes defaultCommsType;

        //- Number of polling cycles in processor updates
        static int nPollProcInterfaces;

        //- Optional maximum message size (bytes)
        static int maxCommsSize;

        //- MPI buffer-size (bytes)
        static const int mpiBufferSize;


    // Standard Communicators

        //- Default world communicator (all processors).
        //- May differ from globalComm if local worlds are in use
        static label worldComm;

        //- Debugging: warn for use of any communicator differing from warnComm
        static label warnComm;

        //- Communicator for all processors, irrespective of any local worlds
        static constexpr label globalComm = 0;

        //- A communicator within the current rank only
        static constexpr label selfComm = 1;

        //- Number of currently defined communicators
        static label nComms() noexcept { return parentComm_.size(); }

        //- True if communicator appears to be user-allocated
        static bool isUserComm(label communicator) noexcept
        {
            return (communicator > worldComm && communicator > selfComm);
        }


    // Constructors

        //- Construct for given communication type
        explicit UPstream(const commsTypes commsType)
        :
            commsType_(commsType)
        {}


    // Member Functions

        //- Allocate a new communicator with subRanks of parent communicator
        static label allocateCommunicator
        (
            //! The parent communicator
            const label parent,

            //! The sub-ranks of parent to use (ignore negative values)
            const labelUList& subRanks,

            //! Call allocatePstreamCommunicator
            const bool doPstream = true
        );

        //- Free a previously allocated communicator.
        //  Ignores placeholder (negative) communicators.
        static void freeCommunicator
        (
            const label communicator,
            const bool doPstream = true
        );

        //- Free all communicators
        static void freeCommunicators(const bool doPstream);


        //- Wrapper class for allocating/freeing communicators
        class communicator
        {
            label comm_;

        public:

            //- No copy construct
            communicator(const communicator&) = delete;

            //- No copy assignment
            void operator=(const communicator&) = delete;

            //- Default construct (a placeholder communicator)
            communicator() noexcept : comm_(-1) {}

            //- Move construct, takes ownership
            communicator(communicator&& c) : comm_(c.comm_) { c.comm_ = -1; }

            //- Allocate a communicator based on given parent
            communicator
            (
                //! The parent communicator
                const label parent,

                //! The sub-ranks of parent to use (ignore negative values)
                const labelUList& subRanks,

                //! Call allocatePstreamCommunicator
                const bool doPstream = true
            )
            :
                comm_(allocateCommunicator(parent, subRanks, doPstream))
            {}

            //- Free allocated communicator and group
            ~communicator() { freeCommunicator(comm_); }

            //- True if communicator is non-negative (ie, was allocated)
            bool good() const noexcept { return (comm_ >= 0); }

            //- The communicator label
            label comm() const noexcept { return comm_; }

            //- Free allocated communicator and group
            void reset() { freeCommunicator(comm_); comm_ = -1; }

            //- Allocate with subRanks of parent communicator
            void reset(label parent, const labelUList& subRanks)
            {
                freeCommunicator(comm_);
                comm_ = allocateCommunicator(parent, subRanks);
            }

            //- Take ownership, free allocated communicator and group.
            void reset(communicator&& c)
            {
                if (comm_ != c.comm_) freeCommunicator(comm_);
                comm_ = c.comm_;
                c.comm_ = -1;
            }

            //- Move assignment, takes ownership
            void operator=(communicator&& c) { reset(std::move(c)); }

            //- Cast to label - the same as comm()
            operator label() const noexcept { return comm_; }
        };


        //- Return physical processor number (i.e. processor number in
        //- worldComm) given communicator and processor
        static int baseProcNo(label comm, int procID);

        //- Return processor number in communicator (given physical processor
        //- number) (= reverse of baseProcNo)
        static label procNo(const label comm, const int baseProcID);

        //- Return processor number in communicator (given processor number
        //- and communicator)
        static label procNo
        (
            const label comm,
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

            //- Number of outstanding requests
            static label nRequests() noexcept;

            //- Truncate outstanding requests to given length
            static void resetRequests(const label n);

            //- Wait until all requests (from start onwards) have finished.
            //  A no-op if parRun() == false
            static void waitRequests(const label start = 0);

            //- Wait until request i has finished.
            //  A no-op if parRun() == false
            //  or for placeholder (negative) request indices
            static void waitRequest(const label i);

            //- Non-blocking comms: has request i finished?
            //  A no-op and returns true if parRun() == false
            //  or for placeholder (negative) request indices
            static bool finishedRequest(const label i);

            static int allocateTag(const char* const msg = nullptr);
            static void freeTag(const int tag, const char* const msg = nullptr);


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

        //- Number of ranks in parallel run (for given communicator)
        //- is 1 for serial run
        static label nProcs(const label communicator = worldComm)
        {
            return procIDs_[communicator].size();
        }

        //- Process index of the master (always 0)
        static constexpr int masterNo() noexcept
        {
            return 0;
        }

        //- Am I the master rank
        static bool master(const label communicator = worldComm)
        {
            return myProcNo_[communicator] == masterNo();
        }

        //- Is this process a sub-rank on the communicator
        static bool is_subrank(const label communicator = worldComm)
        {
            return myProcNo_[communicator] > masterNo();
        }

        //- Number of this process (starting from masterNo() = 0)
        static int myProcNo(const label communicator = worldComm)
        {
            return myProcNo_[communicator];
        }

        //- The parent communicator
        static label parent(const label communicator)
        {
            return parentComm_(communicator);
        }

        //- Process IDs within a given communicator
        static List<int>& procID(const label communicator)
        {
            return procIDs_[communicator];
        }


        // Worlds

            //- All worlds
            static const wordList& allWorlds() noexcept
            {
                return allWorlds_;
            }

            //- The indices into allWorlds for all processes
            static const labelList& worldIDs() noexcept
            {
                return worldIDs_;
            }

            //- My worldID
            static label myWorldID()
            {
                return worldIDs_[myProcNo(globalComm)];
            }

            //- My world
            static const word& myWorld()
            {
                return allWorlds_[worldIDs_[myProcNo(globalComm)]];
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

        //- Communication schedule for linear/tree all-to-master (proc 0).
        //- Chooses based on the value of UPstream::nProcsSimpleSum
        static const List<commsStruct>& whichCommunication
        (
            const label communicator = worldComm
        )
        {
            return
            (
                nProcs(communicator) < nProcsSimpleSum
              ? linearCommunication_[communicator]
              : treeCommunication_[communicator]
            );
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
        //  \return the previous value
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

        //- Exchange integer data with all processors (in the communicator).
        //  \c sendData[proci] is the value to send to proci.
        //  After return recvData contains the data from the other processors.
        static void allToAll
        (
            const UList<int32_t>& sendData,
            UList<int32_t>& recvData,
            const label communicator = worldComm
        );

        //- Exchange integer data with all processors (in the communicator).
        //  \c sendData[proci] is the value to send to proci.
        //  After return recvData contains the data from the other processors.
        static void allToAll
        (
            const UList<int64_t>& sendData,
            UList<int64_t>& recvData,
            const label communicator = worldComm
        );


    // Low-level gather/scatter routines

        #undef  Pstream_CommonRoutines
        #define Pstream_CommonRoutines(Native)                                \
                                                                              \
        /*! \brief Receive identically-sized \c Native data from all ranks */ \
        static void mpiGather                                                 \
        (                                                                     \
            const Native* sendData,                                           \
            int sendCount,                                                    \
            char* recvData,                                                   \
            int recvCount,                                                    \
            const label communicator = worldComm                              \
        );                                                                    \
                                                                              \
        /*! \brief Send identically-sized \c Native data to all ranks */      \
        static void mpiScatter                                                \
        (                                                                     \
            const Native* sendData,                                           \
            int sendCount,                                                    \
            char* recvData,                                                   \
            int recvCount,                                                    \
            const label communicator = worldComm                              \
        );                                                                    \

        Pstream_CommonRoutines(char);


        #undef  Pstream_CommonRoutines
        #define Pstream_CommonRoutines(Native)                                \
                                                                              \
        /*! \brief Receive variable length \c Native data from all ranks */   \
        static void gather                                                    \
        (                                                                     \
            const Native* sendData,                                           \
            int sendCount,  /*!< Ignored on master if recvCount[0] == 0 */    \
            Native* recvData,               /*!< Ignored on non-root rank */  \
            const UList<int>& recvCounts,   /*!< Ignored on non-root rank */  \
            const UList<int>& recvOffsets,  /*!< Ignored on non-root rank */  \
            const label communicator = worldComm                              \
        );                                                                    \
                                                                              \
        /*! \brief Send variable length \c Native data to all ranks */        \
        static void scatter                                                   \
        (                                                                     \
            const Native* sendData,         /*!< Ignored on non-root rank */  \
            const UList<int>& sendCounts,   /*!< Ignored on non-root rank */  \
            const UList<int>& sendOffsets,  /*!< Ignored on non-root rank */  \
            Native* recvData,                                                 \
            int recvCount,                                                    \
            const label communicator = worldComm                              \
        );

        Pstream_CommonRoutines(char);
        Pstream_CommonRoutines(int32_t);
        Pstream_CommonRoutines(int64_t);
        Pstream_CommonRoutines(uint32_t);
        Pstream_CommonRoutines(uint64_t);
        Pstream_CommonRoutines(float);
        Pstream_CommonRoutines(double);

        #undef Pstream_CommonRoutines


    // Gather single, contiguous value(s)

        //- Gather individual values into list locations.
        //  On master list length == nProcs, otherwise zero length.
        //  If called in non-parallel mode,
        //  the returned list length is 1 with localValue.
        template<class T>
        static List<T> listGatherValues
        (
            const T& localValue,
            const label communicator = worldComm
        );

        //- Scatter individual values from list locations.
        //  On master input list length == nProcs, ignored on other procs.
        //  If called in non-parallel mode,
        //  returns the first list element (or zero).
        template<class T>
        static T listScatterValues
        (
            const UList<T>& allValues,
            const label communicator = worldComm
        );


    // Broadcast Functions

        //- Broadcast buffer contents to all processes in communicator.
        //- The sizes must match on all processes.
        //  \return True on success
        static bool broadcast
        (
            char* buf,
            const std::streamsize bufSize,
            const label communicator = worldComm,
            const int rootProcNo = masterNo()
        );


    // Logical reductions

        //- Logical (and) reduction (cf. MPI AllReduce)
        static void reduceAnd
        (
            bool& value,
            const label communicator = worldComm
        );

        //- Logical (or) reduction (cf. MPI AllReduce)
        static void reduceOr
        (
            bool& value,
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
UList<UPstream::commsStruct>::operator[](const label procID);

template<>
const UPstream::commsStruct&
UList<UPstream::commsStruct>::operator[](const label procID) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "UPstreamTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
