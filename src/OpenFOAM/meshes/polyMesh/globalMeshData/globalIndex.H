/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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
    globalIndexI.H
    globalIndex.C
    globalIndexTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef globalIndex_H
#define globalIndex_H

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


    // Private Data

        //- Start of proci. Size is nProcs()+1. (so like CompactListList)
        labelList offsets_;


public:

    // Constructors

        //- Default construct
        globalIndex() = default;

        //- Construct from local max size.
        //  Does communication with default communicator and message tag.
        inline explicit globalIndex(const label localSize);

        //- Construct from local max size.
        //  Does communication with given communicator and message tag
        inline globalIndex
        (
            const label localSize,
            const int tag,          //!< message tag
            const label comm,       //!< communicator
            const bool parallel     //!< use parallel comms
        );

        //- Copy construct from list of labels
        inline explicit globalIndex(const labelUList& offsets);

        //- Move construct from list of labels
        inline explicit globalIndex(labelList&& offsets);

        //- Construct from Istream
        explicit globalIndex(Istream& is);


    // Member Functions

        //- Check for default constructed or global sum == 0
        inline bool empty() const;

        //- Global sum of localSizes
        inline label size() const;

        //- The local sizes
        labelList sizes() const;

        //- Global max of localSizes
        inline label maxSize() const;

        //- Const-access to the offsets
        inline const labelList& offsets() const noexcept;


    // Edit

        //- Write-access to the offsets, for changing after construction
        inline labelList& offsets() noexcept;

        //- Reset from local size.
        //  Does communication with default communicator and message tag.
        inline void reset(const label localSize);

        //- Reset from local size.
        //  Does communication with given communicator and message tag
        void reset
        (
            const label localSize,
            const int tag,          //!< message tag
            const label comm,       //!< communicator
            const bool parallel     //!< use parallel comms
        );


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

            //- Is on local processor
            inline bool isLocal(const label i) const;

            //- From local to global index
            inline label toGlobal(const label i) const;

            //- From local to global index
            inline labelList toGlobal(const labelUList& labels) const;

            //- From local to global index (inplace)
            inline void inplaceToGlobal(labelList& labels) const;

            //- From global to local on current processor.
            //  FatalError if not on local processor.
            inline label toLocal(const label i) const;


        // Global (off-processor) queries

            //- Start of proci data
            inline label offset(const label proci) const;

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
                labelList& labels
            ) const;


            //- From global to local on proci
            inline label toLocal(const label proci, const label i) const;

            //- Which processor does global come from? Binary search.
            inline label whichProcID(const label i) const;


        // Other

            //- Collect data in processor order on master (== procIDs[0]).
            //  Offsets needed on master only.
            template<class Container, class Type>
            static void gather
            (
                const labelUList& offsets,
                const label comm,
                const Container& procIDs,
                const UList<Type>& fld,
                List<Type>& allFld,
                const int tag = UPstream::msgType(),
                const Pstream::commsTypes = Pstream::commsTypes::nonBlocking
            );

            //- Collect indirect data in processor order on master
            //  Offsets needed on master only.
            template<class Type, class Addr>
            static void gather
            (
                const labelUList& offsets,
                const label comm,
                const UList<int>& procIDs,
                const IndirectListBase<Type, Addr>& fld,
                List<Type>& allFld,
                const int tag = UPstream::msgType(),
                const Pstream::commsTypes = Pstream::commsTypes::scheduled
            );

            //- Collect data in processor order on master (== procIDs[0]).
            //  Offsets needed on master only.
            template<class Container, class Type>
            void gather
            (
                const label comm,
                const Container& procIDs,
                const UList<Type>& fld,
                List<Type>& allFld,
                const int tag = UPstream::msgType(),
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::nonBlocking
            ) const
            {
                gather(offsets_, comm, procIDs, fld, allFld, tag, commsType);
            }

            //- Collect data in processor order on master.
            //  Does communication with default communicator and message tag.
            template<class Type>
            void gather
            (
                const UList<Type>& fld,
                List<Type>& allFld,
                const int tag = UPstream::msgType(),
                const Pstream::commsTypes = Pstream::commsTypes::nonBlocking
            ) const;

            //- Collect data indirectly in processor order on master.
            //  Does communication with default communicator and message tag.
            template<class Type, class Addr>
            void gather
            (
                const IndirectListBase<Type, Addr>& fld,
                List<Type>& allFld,
                const int tag = UPstream::msgType(),
                const Pstream::commsTypes = Pstream::commsTypes::scheduled
            ) const;

            //- Collect data in processor order on master.
            //  Does communication with default communicator and message tag.
            template<class Type>
            static void gatherOp
            (
                const UList<Type>& fld,
                List<Type>& allFld,
                const int tag = UPstream::msgType(),
                const Pstream::commsTypes = Pstream::commsTypes::nonBlocking
            );


            //- Inplace collect in processor order on master (== procIDs[0]).
            //- Needs offsets only on master.
            template<class Container, class Type>
            static void gather
            (
                const labelUList& offsets,
                const label comm,
                const Container& procIDs,
                List<Type>& fld,
                const int tag = UPstream::msgType(),
                const Pstream::commsTypes = Pstream::commsTypes::nonBlocking
            );

            //- Inplace collect in processor order on master (== procIDs[0]).
            //- Needs offsets only on master.
            template<class Container, class Type>
            void gather
            (
                const label comm,
                const Container& procIDs,
                List<Type>& fld,
                const int tag = UPstream::msgType(),
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::nonBlocking
            ) const
            {
                gather(offsets_, comm, procIDs, fld, tag, commsType);
            }

            //- Inplace collect data in processor order on master
            //  Does communication with default communicator and message tag.
            //  After the gather, the field is zero-sized on the slaves.
            template<class Type>
            void gather
            (
                List<Type>& fld,
                const int tag = UPstream::msgType(),
                const Pstream::commsTypes = Pstream::commsTypes::nonBlocking
            ) const;

            //- Inplace collect data in processor order on master
            //  Does communication with default communicator and message tag.
            //  After the gather, the field is zero-sized on the slaves.
            template<class Type>
            static void gatherOp
            (
                List<Type>& fld,
                const int tag = UPstream::msgType(),
                const Pstream::commsTypes = Pstream::commsTypes::nonBlocking
            );


            //- Distribute data in processor order. Requires fld to be sized!
            template<class Container, class Type>
            static void scatter
            (
                const labelUList& offsets,
                const label comm,
                const Container& procIDs,
                const UList<Type>& allFld,
                UList<Type>& fld,
                const int tag = UPstream::msgType(),
                const Pstream::commsTypes = Pstream::commsTypes::nonBlocking
            );

            //- Distribute data in processor order. Requires fld to be sized!
            template<class Container, class Type>
            void scatter
            (
                const label comm,
                const Container& procIDs,
                const UList<Type>& allFld,
                UList<Type>& fld,
                const int tag = UPstream::msgType(),
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::nonBlocking
            ) const
            {
                scatter(offsets_, comm, procIDs, allFld, fld, tag, commsType);
            }

            //- Distribute data in processor order. Requires fld to be sized!
            //  Does communication with default communicator and message tag.
            template<class Type>
            void scatter
            (
                const UList<Type>& allFld,
                UList<Type>& fld,
                const int tag = UPstream::msgType(),
                const Pstream::commsTypes = Pstream::commsTypes::nonBlocking
            ) const;

            //- Get (potentially remote) data. Elements required given as
            //  global indices
            template<class Type, class CombineOp>
            void get
            (
                List<Type>& allFld,
                const labelUList& globalIds,
                const CombineOp& cop,
                const label comm = Pstream::worldComm,
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
