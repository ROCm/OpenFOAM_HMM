/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2018 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations

class globalIndex;
class labelRange;

Istream& operator>>(Istream& is, globalIndex& gi);
Ostream& operator<<(Ostream& os, const globalIndex& gi);


/*---------------------------------------------------------------------------*\
                           Class globalIndex Declaration
\*---------------------------------------------------------------------------*/

class globalIndex
{
    // Private data

        //- Start of proci. Size is nProcs()+1. (so like CompactListList)
        labelList offsets_;


public:

    // Constructors

        //- Construct null
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

        //- Check for null constructed or global sum == 0
        inline bool empty() const;

        //- Const-access to the offsets
        inline const labelList& offsets() const;


    // Edit

        //- Write-access to the offsets, for changing after construction
        inline labelList& offsets();

        //- Reset from local size.
        //  Does communication with default communicator and message tag.
        void reset(const label localSize);

        //- Reset from local size.
        //  Does communication with given communicator and message tag
        void reset
        (
            const label localSize,
            const int tag,          //!< message tag
            const label comm,       //!< communicator
            const bool parallel     //!< use parallel comms
        );


        // Queries relating to my processor (using world communicator)

            //- My local start
            inline label localStart() const;

            //- My local size
            inline label localSize() const;

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


        // Global queries

            //- Global sum of localSizes
            inline label size() const;

            //- Start of proci data
            inline label offset(const label proci) const;

            //- Start of proci data
            inline label localStart(const label proci) const;

            //- Size of proci data
            inline label localSize(const label proci) const;

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
            //  Needs offsets only on master.
            template<class Type>
            static void gather
            (
                const labelUList& offsets,
                const label comm,
                const labelList& procIDs,
                const UList<Type>& fld,
                List<Type>& allFld,
                const int tag = UPstream::msgType(),
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::nonBlocking
            );

            //- Collect data in processor order on master (== procIDs[0]).
            //  Needs offsets only on master.
            template<class Type>
            void gather
            (
                const label comm,
                const labelList& procIDs,
                const UList<Type>& fld,
                List<Type>& allFld,
                const int tag = UPstream::msgType(),
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::nonBlocking
            ) const
            {
                gather(offsets_, comm, procIDs, fld, allFld, tag, commsType);
            }

            //- Inplace collect data in processor order on master
            //  (== procIDs[0]). Needs offsets only on master.
            template<class Type>
            static void gather
            (
                const labelUList& offsets,
                const label comm,
                const labelList& procIDs,
                List<Type>& fld,
                const int tag = UPstream::msgType(),
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::nonBlocking
            );

            //- Inplace collect data in processor order on master
            //  (== procIDs[0]). Needs offsets only on master.
            template<class Type>
            void gather
            (
                const label comm,
                const labelList& procIDs,
                List<Type>& fld,
                const int tag = UPstream::msgType(),
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::nonBlocking
            ) const
            {
                gather(offsets_, comm, procIDs, fld, tag, commsType);
            }

            //- Distribute data in processor order. Requires fld to be sized!
            template<class Type>
            static void scatter
            (
                const labelUList& offsets,
                const label comm,
                const labelList& procIDs,
                const UList<Type>& allFld,
                UList<Type>& fld,
                const int tag = UPstream::msgType(),
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::nonBlocking
            );

            //- Distribute data in processor order. Requires fld to be sized!
            template<class Type>
            void scatter
            (
                const label comm,
                const labelList& procIDs,
                const UList<Type>& allFld,
                UList<Type>& fld,
                const int tag = UPstream::msgType(),
                const Pstream::commsTypes commsType =
                    Pstream::commsTypes::nonBlocking
            ) const
            {
                scatter(offsets_, comm, procIDs, allFld, fld, tag, commsType);
            }


    // IOstream Operators

        friend Istream& operator>>(Istream& is, globalIndex& gi);
        friend Ostream& operator<<(Ostream& os, const globalIndex& gi);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "globalIndexTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "globalIndexI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
