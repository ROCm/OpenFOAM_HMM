/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2022 OpenCFD Ltd.
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
    Foam::globalIndex

Description
    Calculates a unique integer (label so might not have enough room - 2G max)
    for processor + local index. E.g.

    globalIndex globalFaces(mesh.nFaces());
    label globalFacei = globalFaces.toGlobal(facei);

SourceFiles
    globalIndex.C
    globalIndexI.H
    globalIndexTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_globalIndex_H
#define Foam_globalIndex_H

#include "Pstream.H"
#include "CompactListList.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class globalIndex;
class labelRange;

Istream& operator>>(Istream& is, globalIndex& gi);
Ostream& operator<<(Ostream& os, const globalIndex& gi);


/*---------------------------------------------------------------------------*\
                           Class globalIndex Declaration
\*---------------------------------------------------------------------------*/

class globalIndex
{
    // Private Data

        //- Start of proci. Size is nProcs()+1. (so like CompactListList)
        labelList offsets_;


    // Private Member Functions

        //- Sort and bin. validBins contains bins with non-zero size.
        static void bin
        (
            const labelUList& offsets,
            const labelUList& globalIds,
            labelList& order,
            CompactListList<label>& sortedElems,
            DynamicList<label>& validBins
        );

        //- Report overflow at specified index
        static void reportOverflowAndExit
        (
            const label idx,
            const labelUList& localLens = labelUList::null()
        );

public:

    // Public Data Types

        //- Disambiguation tag (list construction dispatch)
        enum accessType : char { OFFSETS, SIZES };

        //- Dispatch tag
        struct gatherOnly{};

        //- Dispatch tag
        struct gatherNone{};


    // Constructors

        //- Default construct
        globalIndex() = default;

        //- Copy construct from a list of offsets.
        //- No communication required
        inline explicit globalIndex(const labelUList& listOffsets);

        //- Move construct from a list of offsets.
        //- No communication required
        inline explicit globalIndex(labelList&& listOffsets);

        //- Copy construct from a list of offsets or sizes.
        //- No communication required
        inline globalIndex
        (
            const labelUList& offsetsOrSizes,
            enum accessType accType
        );

        //- Construct from local size.
        //  Communication with default communicator and message tag.
        inline explicit globalIndex(const label localSize);

        //- Construct by gathering local sizes without rescattering.
        //- This 'one-sided' globalIndex will be empty on non-master processes.
        //
        //  \note can be used when Pstream::parRun() is false.
        inline globalIndex
        (
            const label localSize,
            const globalIndex::gatherOnly,
            const label comm = UPstream::worldComm  //!< communicator
        );

        //- Construct with a single size entry.
        //- No communication required
        inline globalIndex
        (
            const label localSize,
            const globalIndex::gatherNone,
            const label comm = -1   //!< dummy communicator
        );

        //- Construct from local size.
        //  Communication with given communicator and message tag,
        //  unless parallel == false
        inline globalIndex
        (
            const label localSize,
            const int tag,          //!< message tag
            const label comm,       //!< communicator
            const bool parallel     //!< use parallel comms
        );

        //- Construct from Istream.
        //- No communication required
        explicit globalIndex(Istream& is);


    // Member Functions

        //- Check for default constructed or global sum == 0
        inline bool empty() const;

        //- Global sum of localSizes. Same as totalSize()
        inline label size() const;

        //- Global sum of localSizes.
        inline label totalSize() const;

        //- The local sizes. Same as localSizes()
        inline labelList sizes() const;

        //- The local starts
        inline const labelUList localStarts() const;

        //- The local sizes
        labelList localSizes() const;

        //- Global max of localSizes
        inline label maxSize() const;

        //- Const-access to the offsets
        inline const labelList& offsets() const noexcept;


    // Dimensions

        //- The number of processors covered by the offsets
        inline label nProcs() const noexcept;

        //- Range of process indices for all addressed offsets (processes)
        inline labelRange allProcs() const noexcept;

        //- Range of process indices for addressed sub-offsets (processes)
        inline labelRange subProcs() const noexcept;


    // Edit

        //- Write-access to the offsets, for changing after construction
        inline labelList& offsets() noexcept;

        //- Reset from local size.
        //  Does communication with default communicator and message tag.
        void reset(const label localSize);

        //- Reset by gathering local sizes without rescattering.
        //- This 'one-sided' globalIndex will be empty on non-master processes.
        //
        //  \note can be used when Pstream::parRun() is false.
        void reset
        (
            const label localSize,
            const globalIndex::gatherOnly,
            const label comm = UPstream::worldComm  //!< communicator
        );

        //- Reset from local size.
        //  Does communication with given communicator and message tag,
        //  unless parallel == false
        void reset
        (
            const label localSize,
            const int tag,          //!< message tag
            const label comm,       //!< communicator
            const bool parallel     //!< use parallel comms
        );

        //- Reset from list of local sizes,
        //- with optional check for label overflow.
        //- No communication required
        void reset
        (
            const labelUList& sizes,
            const bool checkOverflow = false
        );

        //- Alter local size for given processor
        void setLocalSize(const label proci, const label len);


    // Queries

        // Queries relating to my processor (using world communicator)

            //- My local start
            inline label localStart() const;

            //- My local size
            inline label localSize() const;

            //- The max of localSizes, excluding current processor
            inline label maxNonLocalSize() const;

            //- Return start/size range of local processor data
            inline labelRange range() const;

            //- Return start/size ranges for all data
            List<labelRange> ranges() const;

            //- Is on local processor
            inline bool isLocal(const label i) const;

            //- From local to global index
            inline label toGlobal(const label i) const;

            //- From local to global index
            inline labelList toGlobal(const labelUList& labels) const;

            //- From local to global index (inplace)
            inline void inplaceToGlobal(labelUList& labels) const;

            //- From global to local on current processor.
            //  FatalError if not on local processor.
            inline label toLocal(const label i) const;


        // Global (off-processor) queries

            //- Start of proci data
            inline label localStart(const label proci) const;

            //- Size of proci data
            inline label localSize(const label proci) const;

            //- The max of localSizes, excluding the specified processor
            label maxNonLocalSize(const label proci) const;

            //- Return start/size range of proci data
            inline labelRange range(const label proci) const;

            //- Is on processor proci
            inline bool isLocal(const label proci, const label i) const;

            //- From local to global on proci
            inline label toGlobal(const label proci, const label i) const;

            //- From local to global on proci
            inline labelList toGlobal
            (
                const label proci,
                const labelUList& labels
            ) const;

            //- From local to global index on proci (inplace)
            inline void inplaceToGlobal
            (
                const label proci,
                labelUList& labels
            ) const;


            //- From global to local on proci
            inline label toLocal(const label proci, const label i) const;

            //- Which processor does global come from? Binary search.
            inline label whichProcID(const label i) const;


    // Housekeeping

            //- Same as localStart
            label offset(const label proci) const
            {
                return localStart(proci);
            }


    // Iteration

        //- Forward input iterator with const access
        class const_iterator
        {
            //- The parent class for which this is an iterator
            const globalIndex* parent_;

            //- The index into the parent
            label proci_;


        public:

        // Constructors

            //- Construct from globalIndex list at given index
            explicit const_iterator
            (
                const globalIndex* globalIdx,
                const label proci = 0
            ) noexcept;


        // Member Operators

            //- The associated local proc
            inline label proci() const noexcept;

            //- The local start
            inline label start() const;

            //- The local size
            inline label size() const;

            //- The local range
            inline labelRange range() const;

            //- The local range
            inline labelRange operator*() const;

            inline const_iterator& operator++();
            inline const_iterator operator++(int);

            inline const_iterator& operator--();
            inline const_iterator operator--(int);

            inline bool operator==(const const_iterator& iter) const noexcept;
            inline bool operator!=(const const_iterator& iter) const noexcept;
        };


        //- A const_iterator set to the beginning
        inline const_iterator cbegin() const noexcept;

        //- A const_iterator set to beyond the end
        inline const const_iterator cend() const noexcept;

        //- A const_iterator set to the beginning
        inline const_iterator begin() const noexcept;

        //- A const_iterator set to beyond the end
        inline const const_iterator end() const noexcept;


    // Helper Functions

        //- Calculate offsets from a list of local sizes,
        //- with optional check for label overflow
        static labelList calcOffsets
        (
            const labelUList& localLens,
            const bool checkOverflow = false
        );

        //- Calculate offsets from list of lists,
        //- with optional check for label overflow
        template<class SubListType>
        static labelList calcListOffsets
        (
            const List<SubListType>& lists,
            const bool checkOverflow = false
        );

        //- Calculate ranges (offset/size) from a list of local sizes,
        //- with optional check for label overflow
        static List<labelRange> calcRanges
        (
            const labelUList& localLens,
            const bool checkOverflow = false
        );


    // Other

        // Gather (helpers)

            //- Collect data in processor order on master (== procIDs[0]).
            //  \note offsets needed on master only.
            template<class ProcIDsContainer, class Type>
            static void gather
            (
                const labelUList& offsets,
                const label comm,       //!< communicator
                const ProcIDsContainer& procIDs,
                const UList<Type>& fld,
                List<Type>& allFld,
                const int tag = UPstream::msgType(),
                const Pstream::commsTypes = Pstream::commsTypes::nonBlocking
            );

            //- Collect indirect data in processor order on master
            //  \note offsets needed on master only.
            template<class Type, class Addr>
            static void gather
            (
                const labelUList& offsets,
                const label comm,       //!< communicator
                const UList<int>& procIDs,
                const IndirectListBase<Type, Addr>& fld,
                List<Type>& allFld,
                const int tag = UPstream::msgType(),
                const Pstream::commsTypes = Pstream::commsTypes::scheduled
            );

            //- Inplace collect in processor order on master (== procIDs[0]).
            //  \note offsets needed on master only.
            template<class ProcIDsContainer, class Type>
            static void gather
            (
                const labelUList& offsets,
                const label comm,       //!< communicator
                const ProcIDsContainer& procIDs,
                List<Type>& fld,
                const int tag = UPstream::msgType(),
                const Pstream::commsTypes = Pstream::commsTypes::nonBlocking
            );

            //- Collect data in processor order on master (== procIDs[0]).
            //  \note the globalIndex offsets needed on master only.
            template<class ProcIDsContainer, class Type>
            void gather
            (
                const label comm,       //!< communicator
                const ProcIDsContainer& procIDs,
                const UList<Type>& fld,
                List<Type>& allFld,
                const int tag = UPstream::msgType(),
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::nonBlocking
            ) const
            {
                gather(offsets_, comm, procIDs, fld, allFld, tag, commsType);
            }

            //- Inplace collect in processor order on master (== procIDs[0]).
            //  \note the globalIndex offsets needed on master only.
            template<class ProcIDsContainer, class Type>
            void gather
            (
                const label comm,       //!< communicator
                const ProcIDsContainer& procIDs,
                List<Type>& fld,
                const int tag = UPstream::msgType(),
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::nonBlocking
            ) const
            {
                gather(offsets_, comm, procIDs, fld, tag, commsType);
            }


        // Gather

            //- Collect data in processor order on master
            //- (in serial: performs a simple copy).
            //  Communication with default/specified communicator, message tag.
            template<class Type>
            void gather
            (
                const UList<Type>& sendData,
                List<Type>& allData,
                const int tag = UPstream::msgType(),
                const Pstream::commsTypes = Pstream::commsTypes::nonBlocking,
                const label comm = UPstream::worldComm  //!< communicator
            ) const;

            //- Collect data indirectly in processor order on master
            //- (in serial: performs a simple copy).
            //  Communication with default/specified communicator, message tag.
            template<class Type, class Addr>
            void gather
            (
                const IndirectListBase<Type, Addr>& sendData,
                List<Type>& allData,
                const int tag = UPstream::msgType(),
                const Pstream::commsTypes = Pstream::commsTypes::scheduled,
                const label comm = UPstream::worldComm  //!< communicator
            ) const;

            //- Collect data in processor order on master
            //- (in serial: performs a simple copy).
            //  Communication with default/specified communicator, message tag.
            //
            //  \return output (master), zero-sized on non-master
            template<class Type, class OutputContainer = List<Type>>
            OutputContainer gather
            (
                const UList<Type>& sendData,
                const int tag = UPstream::msgType(),
                const Pstream::commsTypes = Pstream::commsTypes::nonBlocking,
                const label comm = UPstream::worldComm  //!< communicator
            ) const;

            //- Collect data indirectly in processor order on master.
            //  Communication with default/specified communicator, message tag.
            //
            //  \return output (master), zero-sized on non-master
            template<class Type, class Addr, class OutputContainer = List<Type>>
            OutputContainer gather
            (
                const IndirectListBase<Type, Addr>& sendData,
                const int tag = UPstream::msgType(),
                const Pstream::commsTypes = Pstream::commsTypes::scheduled,
                const label comm = UPstream::worldComm  //!< communicator
            ) const;

            //- Inplace collect data in processor order on master
            //- (in serial: a no-op).
            //  Communication with default/specified communicator, message tag.
            //  After the gather, the field is zero-sized on the slaves.
            template<class Type>
            void gatherInplace
            (
                List<Type>& fld,
                const int tag = UPstream::msgType(),
                const Pstream::commsTypes = Pstream::commsTypes::nonBlocking,
                const label comm = UPstream::worldComm  //!< communicator
            ) const;

            //- Collect \em contiguous data using a MPI_Gatherv call
            //- (in serial: performs a simple copy).
            //  Communication with default/specified communicator.
            //  \attention The nProcs for globalIndex and communicator
            //      must match!!
            //
            //  The allData is output (master), zero-sized on non-master
            template<class Type, class OutputContainer = List<Type>>
            void mpiGather
            (
                const UList<Type>& sendData,
                OutputContainer& allData,
                const label comm = UPstream::worldComm  //!< communicator
            ) const;

            //- Collect \em contiguous data using a MPI_Gatherv call
            //- (in serial: performs a simple copy).
            //  Communication with default/specified communicator.
            //  \attention The nProcs for globalIndex and communicator
            //      must match!!
            //
            //  \return output (master), zero-sized on non-master
            template<class Type, class OutputContainer = List<Type>>
            OutputContainer mpiGather
            (
                const UList<Type>& sendData,
                const label comm = UPstream::worldComm  //!< communicator
            ) const;

            //- Inplace collect \em contiguous data using a MPI_Gatherv call
            //- (in serial: a no-op).
            //  Communication with default/specified communicator.
            //  \attention The nProcs for globalIndex and communicator
            //      must match!!
            //
            //  After the gather, the field is zero-sized on non-master.
            template<class Type>
            void mpiGatherInplace
            (
                List<Type>& fld,
                const label comm = UPstream::worldComm  //!< communicator
            ) const;


        // Gather Operations

            //- Collect \em contiguous data using a MPI_Gatherv call
            //- (in serial: performs a simple copy).
            //  Communication with default/specified communicator.
            //
            //  The allData is output (master), zero-sized on non-master
            template<class Type, class OutputContainer = List<Type>>
            static void mpiGatherOp
            (
                const UList<Type>& sendData,
                OutputContainer& allData,
                const label comm = UPstream::worldComm  //!< communicator
            );

            //- Collect \em contiguous data using a MPI_Gatherv call
            //- (in serial: performs a simple copy).
            //  Communication with default/specified communicator.
            //
            //  \return output (master), zero-sized on non-master
            template<class Type, class OutputContainer = List<Type>>
            static OutputContainer mpiGatherOp
            (
                const UList<Type>& sendData,
                const label comm = UPstream::worldComm  //!< communicator
            );

            //- Inplace collect \em contiguous data using a MPI_Gatherv call
            //- (in serial: a no-op).
            //  Communication with default/specified communicator.
            //
            //  After the gather, the field is zero-sized on non-master.
            template<class Type>
            static void mpiGatherInplaceOp
            (
                List<Type>& fld,
                const label comm = UPstream::worldComm  //!< communicator
            );

            //- Collect data in processor order on master
            //- (in serial: performs a simple copy).
            //  Communication with default/specified communicator, message tag.
            //
            //  The allFld is output (master), zero-sized on non-master
            template<class Type>
            static void gatherOp
            (
                const UList<Type>& sendData,
                List<Type>& allData,
                const int tag = UPstream::msgType(),
                const Pstream::commsTypes = Pstream::commsTypes::nonBlocking,
                const label comm = UPstream::worldComm  //!< communicator
            );

            //- Collect data in processor order on master
            //- (in serial: performs a simple copy).
            //  Communication with default/specified communicator, message tag.
            //
            //  The allFld is output (master), zero-sized on non-master
            template<class Type, class Addr>
            static void gatherOp
            (
                const IndirectListBase<Type, Addr>& sendData,
                List<Type>& allData,
                const int tag = UPstream::msgType(),
                const Pstream::commsTypes = Pstream::commsTypes::nonBlocking,
                const label comm = UPstream::worldComm  //!< communicator
            );

            //- Collect and return data in processor order on master
            //- (in serial: performs a simple copy).
            //  Communication with default/specified communicator, message tag.
            //
            //  \return output (master), zero-sized on non-master
            template<class Type, class OutputContainer = List<Type>>
            static OutputContainer gatherOp
            (
                const UList<Type>& sendData,
                const int tag = UPstream::msgType(),
                const Pstream::commsTypes = Pstream::commsTypes::nonBlocking,
                const label comm = UPstream::worldComm  //!< communicator
            );

            //- Collect and return data in processor order on master
            //- (in serial: performs a simple copy).
            //  Communication with default/specified communicator, message tag.
            //
            //  \return output (master), zero-sized on non-master
            template<class Type, class Addr, class OutputContainer = List<Type>>
            static OutputContainer gatherOp
            (
                const IndirectListBase<Type, Addr>& sendData,
                const int tag = UPstream::msgType(),
                const Pstream::commsTypes = Pstream::commsTypes::nonBlocking,
                const label comm = UPstream::worldComm  //!< communicator
            );

            //- Inplace collect data in processor order on master
            //- (in serial: a no-op).
            //  Communication with default/specified communicator, message tag.
            //
            //  After the gather, the field is zero-sized on the slaves.
            template<class Type>
            static void gatherInplaceOp
            (
                List<Type>& fld,
                const int tag = UPstream::msgType(),
                const Pstream::commsTypes = Pstream::commsTypes::nonBlocking,
                const label comm = UPstream::worldComm  //!< communicator
            );


        // Scatter

            //- Distribute data in processor order.
            //  Requires fld to be correctly sized!
            //  Communication with default/specified communicator, message tag.
            //  \note offsets needed on master only.
            template<class ProcIDsContainer, class Type>
            static void scatter
            (
                const labelUList& offsets,
                const label comm,       //!< communicator
                const ProcIDsContainer& procIDs,
                const UList<Type>& allFld,
                UList<Type>& fld,
                const int tag = UPstream::msgType(),
                const Pstream::commsTypes = Pstream::commsTypes::nonBlocking
            );

            //- Distribute data in processor order.
            //  Requires fld to be correctly sized!
            //  Communication with default/specified communicator, message tag.
            //  \note the globalIndex offsets needed on master only.
            template<class ProcIDsContainer, class Type>
            void scatter
            (
                const label comm,       //!< communicator
                const ProcIDsContainer& procIDs,
                const UList<Type>& allFld,
                UList<Type>& fld,
                const int tag = UPstream::msgType(),
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::nonBlocking
            ) const
            {
                scatter(offsets_, comm, procIDs, allFld, fld, tag, commsType);
            }

            //- Distribute data in processor order.
            //  Requires fld to be correctly sized!
            //  Communication with default/specified communicator, message tag.
            //  \note the globalIndex offsets needed on master only.
            template<class Type>
            void scatter
            (
                const UList<Type>& allData,
                UList<Type>& localData,
                const int tag = UPstream::msgType(),
                const Pstream::commsTypes = Pstream::commsTypes::nonBlocking,
                const label comm = UPstream::worldComm  //!< communicator
            ) const;

            //- Distribute data in processor order
            //- (in serial: performs a simple copy).
            //  Communication with default/specified communicator, message tag.
            //  \note the globalIndex offsets needed on master only.
            template<class Type, class OutputContainer = List<Type>>
            OutputContainer scatter
            (
                const UList<Type>& allData,
                const int tag = UPstream::msgType(),
                const Pstream::commsTypes = Pstream::commsTypes::nonBlocking,
                const label comm = UPstream::worldComm  //!< communicator
            ) const;


        // Scatter

            //- Get (potentially remote) data.
            //- Elements required given as global indices
            //  Communication with default/specified communicator, message tag.
            template<class Type, class CombineOp>
            void get
            (
                List<Type>& allFld,
                const labelUList& globalIds,
                const CombineOp& cop,
                const label comm = UPstream::worldComm,  //!< communicator
                const int tag = UPstream::msgType()
            ) const;


    // IOstream Operators

        friend Istream& operator>>(Istream& is, globalIndex& gi);
        friend Ostream& operator<<(Ostream& os, const globalIndex& gi);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "globalIndexI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "globalIndexTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
