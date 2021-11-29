/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "globalIndex.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class SubListType>
Foam::labelList
Foam::globalIndex::calcListOffsets
(
    const List<SubListType>& lists,
    const bool checkOverflow
)
{
    labelList values;

    const label len = lists.size();

    if (len)
    {
        values.resize(len+1);

        label start = 0;
        for (label i = 0; i < len; ++i)
        {
            values[i] = start;
            start += lists[i].size();

            if (checkOverflow && start < values[i])
            {
                FatalErrorInFunction
                    << "Overflow : sum of sizes exceeds labelMax ("
                    << labelMax << ") after index " << i << nl
                    << "Please recompile with larger datatype for label." << nl
                    << exit(FatalError);
            }
        }
        values[len] = start;
    }

    return values;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ProcIDsContainer, class Type>
void Foam::globalIndex::gather
(
    const labelUList& off,
    const label comm,
    const ProcIDsContainer& procIDs,
    const UList<Type>& fld,
    List<Type>& allFld,
    const int tag,
    const Pstream::commsTypes commsType
)
{
    if
    (
        !is_contiguous<Type>::value
     && commsType == Pstream::commsTypes::nonBlocking
    )
    {
        FatalErrorInFunction
            << "Cannot use nonBlocking with non-contiguous data"
            << exit(FatalError);
        // Could also warn and change to scheduled etc...
    }

    if (Pstream::myProcNo(comm) == procIDs[0])
    {
        allFld.resize_nocopy(off.last());

        // Assign my local data - respect offset information
        // so that we can request 0 entries to be copied.
        // Also handle the case where we have a slice of the full
        // list.

        SubList<Type>(allFld, off[1]-off[0], off[0]) =
           SubList<Type>(fld, off[1]-off[0]);

        if
        (
            commsType == Pstream::commsTypes::scheduled
         || commsType == Pstream::commsTypes::blocking
        )
        {
            for (label i = 1; i < procIDs.size(); ++i)
            {
                SubList<Type> procSlot(allFld, off[i+1]-off[i], off[i]);

                if (is_contiguous<Type>::value)
                {
                    IPstream::read
                    (
                        commsType,
                        procIDs[i],
                        procSlot.data_bytes(),
                        procSlot.size_bytes(),
                        tag,
                        comm
                    );
                }
                else
                {
                    IPstream fromProc
                    (
                        commsType,
                        procIDs[i],
                        0,
                        tag,
                        comm
                    );
                    fromProc >> procSlot;
                }
            }
        }
        else
        {
            // nonBlocking && is_contiguous == true (already checked)

            const label startOfRequests = Pstream::nRequests();

            // Set up reads
            for (label i = 1; i < procIDs.size(); ++i)
            {
                SubList<Type> procSlot(allFld, off[i+1]-off[i], off[i]);

                IPstream::read
                (
                    commsType,
                    procIDs[i],
                    procSlot.data_bytes(),
                    procSlot.size_bytes(),
                    tag,
                    comm
                );
            }

            // Wait for all to finish
            Pstream::waitRequests(startOfRequests);
        }
    }
    else
    {
        if
        (
            commsType == Pstream::commsTypes::scheduled
         || commsType == Pstream::commsTypes::blocking
        )
        {
            if (is_contiguous<Type>::value)
            {
                OPstream::write
                (
                    commsType,
                    procIDs[0],
                    fld.cdata_bytes(),
                    fld.size_bytes(),
                    tag,
                    comm
                );
            }
            else
            {
                OPstream toMaster
                (
                    commsType,
                    procIDs[0],
                    0,
                    tag,
                    comm
                );
                toMaster << fld;
            }
        }
        else
        {
            // nonBlocking && is_contiguous == true (already checked)

            const label startOfRequests = Pstream::nRequests();

            // Set up write
            OPstream::write
            (
                commsType,
                procIDs[0],
                fld.cdata_bytes(),
                fld.size_bytes(),
                tag,
                comm
            );

            // Wait for all to finish
            Pstream::waitRequests(startOfRequests);
        }
    }
}


template<class Type, class Addr>
void Foam::globalIndex::gather
(
    const labelUList& off,
    const label comm,
    const UList<int>& procIDs,
    const IndirectListBase<Type, Addr>& fld,
    List<Type>& allFld,
    const int tag,
    const Pstream::commsTypes commsType
)
{
    if (commsType == Pstream::commsTypes::nonBlocking)
    {
        WarningInFunction
            << "Cannot use nonBlocking with indirect list of data"
            << exit(FatalError);
        // Could also warn and change to scheduled etc...
    }

    if (Pstream::myProcNo(comm) == procIDs[0])
    {
        allFld.resize_nocopy(off.last());

        // Assign my local data - respect offset information
        // so that we can request 0 entries to be copied

        SubList<Type> localSlot(allFld, off[1]-off[0], off[0]);
        if (!localSlot.empty())
        {
            localSlot = fld;
        }

        // Already verified commsType != nonBlocking
        for (label i = 1; i < procIDs.size(); ++i)
        {
            SubList<Type> procSlot(allFld, off[i+1]-off[i], off[i]);

            IPstream fromProc
            (
                commsType,
                procIDs[i],
                0,
                tag,
                comm
            );
            fromProc >> procSlot;
        }
    }
    else
    {
        OPstream toMaster
        (
            commsType,
            procIDs[0],
            0,
            tag,
            comm
        );
        toMaster << fld;
    }
}


template<class Type>
void Foam::globalIndex::gather
(
    const UList<Type>& fld,
    List<Type>& allFld,
    const int tag,
    const Pstream::commsTypes commsType,
    const label comm
) const
{
    gather
    (
        comm,
        UPstream::procID(comm),
        fld,
        allFld,
        tag,
        commsType
    );
}


template<class Type, class Addr>
void Foam::globalIndex::gather
(
    const IndirectListBase<Type, Addr>& fld,
    List<Type>& allFld,
    const int tag,
    const Pstream::commsTypes commsType,
    const label comm
) const
{
    gather
    (
        offsets_,
        comm,
        UPstream::procID(comm),
        fld,
        allFld,
        tag,
        commsType
    );
}


template<class ProcIDsContainer, class Type>
void Foam::globalIndex::gather
(
    const labelUList& off,
    const label comm,
    const ProcIDsContainer& procIDs,
    List<Type>& fld,
    const int tag,
    const Pstream::commsTypes commsType
)
{
    List<Type> allFld;

    gather(off, comm, procIDs, fld, allFld, tag, commsType);

    if (Pstream::myProcNo(comm) == procIDs[0])
    {
        fld.transfer(allFld);
    }
}


template<class Type>
void Foam::globalIndex::gather
(
    List<Type>& fld,
    const int tag,
    const Pstream::commsTypes commsType,
    const label comm
) const
{
    List<Type> allFld;

    gather
    (
        comm,
        UPstream::procID(comm),
        fld,
        allFld,
        tag,
        commsType
    );

    if (Pstream::master(comm))
    {
        fld.transfer(allFld);
    }
    else
    {
        fld.clear();
    }
}


template<class Type, class OutputContainer>
void Foam::globalIndex::mpiGather
(
    const UList<Type>& sendData,
    OutputContainer& allValues,
    const label comm
) const
{
    if (!is_contiguous<Type>::value)
    {
        FatalErrorInFunction
            << "Cannot be called for non-contiguous data" << nl
            << abort(FatalError);
    }

    const label proci = Pstream::myProcNo(comm);

    const globalIndex& globalAddr = *this;

    // Must be the same as Pstream::nProcs(comm), at least on master!!
    const label nproc = globalAddr.nProcs();

    auto nSendBytes = sendData.size_bytes();

    // Respect local size information so that we can request
    // 0 entries to be sent on master

    if (proci < nproc && !globalAddr.localSize(proci))
    {
        nSendBytes = 0;
    }

    List<int> recvSizes;
    List<int> recvOffsets;

    if (Pstream::master(comm))
    {
        allValues.resize_nocopy(globalAddr.size());

        recvSizes.resize(nproc);
        recvOffsets.resize(nproc+1);

        for (label proci = 0; proci < nproc; ++proci)
        {
            recvSizes[proci] = globalAddr.localSize(proci) * sizeof(Type);
            recvOffsets[proci] = globalAddr.localStart(proci) * sizeof(Type);
        }
        recvOffsets[nproc] = globalAddr.size() * sizeof(Type);
    }
    else
    {
        allValues.clear();
    }

    UPstream::gather
    (
        sendData.cdata_bytes(),
        nSendBytes,
        allValues.data_bytes(),
        recvSizes,
        recvOffsets,
        comm
    );
}


template<class Type, class OutputContainer>
OutputContainer Foam::globalIndex::mpiGather
(
    const UList<Type>& sendData,
    const label comm
) const
{
    OutputContainer allValues;
    mpiGather<Type, OutputContainer>(sendData, allValues, comm);
    return allValues;
}


template<class Type>
void Foam::globalIndex::gatherOp
(
    const UList<Type>& fld,
    List<Type>& allFld,
    const int tag,
    const Pstream::commsTypes commsType
)
{
    globalIndex(fld.size()).gather(fld, allFld, tag, commsType);
}


template<class Type>
void Foam::globalIndex::gatherOp
(
    List<Type>& fld,
    const int tag,
    const Pstream::commsTypes commsType
)
{
    globalIndex(fld.size()).gather(fld, tag, commsType);
}


template<class ProcIDsContainer, class Type>
void Foam::globalIndex::scatter
(
    const labelUList& off,
    const label comm,
    const ProcIDsContainer& procIDs,
    const UList<Type>& allFld,
    UList<Type>& fld,
    const int tag,
    const Pstream::commsTypes commsType
)
{
    if
    (
        !is_contiguous<Type>::value
     && commsType == Pstream::commsTypes::nonBlocking
    )
    {
        FatalErrorInFunction
            << "Cannot use nonBlocking with non-contiguous data"
            << exit(FatalError);
        // Could also warn and change to scheduled etc...
    }

    if (Pstream::myProcNo(comm) == procIDs[0])
    {
        const SubList<Type> localSlot(allFld, off[1]-off[0], off[0]);

        if (!localSlot.empty())
        {
            fld.deepCopy(localSlot);
        }

        if
        (
            commsType == Pstream::commsTypes::scheduled
         || commsType == Pstream::commsTypes::blocking
        )
        {
            for (label i = 1; i < procIDs.size(); ++i)
            {
                const SubList<Type> procSlot(allFld, off[i+1]-off[i], off[i]);

                if (is_contiguous<Type>::value)
                {
                    OPstream::write
                    (
                        commsType,
                        procIDs[i],
                        procSlot.cdata_bytes(),
                        procSlot.size_bytes(),
                        tag,
                        comm
                    );
                }
                else
                {
                    OPstream toProc
                    (
                        commsType,
                        procIDs[i],
                        0,
                        tag,
                        comm
                    );
                    toProc << procSlot;
                }
            }
        }
        else
        {
            // nonBlocking && is_contiguous == true (already checked)

            const label startOfRequests = Pstream::nRequests();

            // Set up writes
            for (label i = 1; i < procIDs.size(); ++i)
            {
                const SubList<Type> procSlot(allFld, off[i+1]-off[i], off[i]);

                OPstream::write
                (
                    commsType,
                    procIDs[i],
                    procSlot.cdata_bytes(),
                    procSlot.size_bytes(),
                    tag,
                    comm
                );
            }

            // Wait for all to finish
            Pstream::waitRequests(startOfRequests);
        }
    }
    else
    {
        if
        (
            commsType == Pstream::commsTypes::scheduled
         || commsType == Pstream::commsTypes::blocking
        )
        {
            if (is_contiguous<Type>::value)
            {
                IPstream::read
                (
                    commsType,
                    procIDs[0],
                    fld.data_bytes(),
                    fld.size_bytes(),
                    tag,
                    comm
                );
            }
            else
            {
                IPstream fromMaster
                (
                    commsType,
                    procIDs[0],
                    0,
                    tag,
                    comm
                );
                fromMaster >> fld;
            }
        }
        else
        {
            // nonBlocking && is_contiguous == true (already checked)

            const label startOfRequests = Pstream::nRequests();

            // Set up read
            IPstream::read
            (
                commsType,
                procIDs[0],
                fld.data_bytes(),
                fld.size_bytes(),
                tag,
                comm
            );

            // Wait for all to finish
            Pstream::waitRequests(startOfRequests);
        }
    }
}


template<class Type>
void Foam::globalIndex::scatter
(
    const UList<Type>& allFld,
    UList<Type>& fld,
    const int tag,
    const Pstream::commsTypes commsType,
    const label comm
) const
{
    scatter
    (
        offsets_,
        comm,
        UPstream::procID(comm),
        allFld,
        fld,
        tag,
        commsType
    );
}


template<class Type, class CombineOp>
void Foam::globalIndex::get
(
    List<Type>& allFld,
    const labelUList& globalIds,
    const CombineOp& cop,
    const label comm,
    const int tag
) const
{
    allFld.resize_nocopy(globalIds.size());
    if (globalIds.size())
    {
        // Sort according to processor
        labelList order;
        CompactListList<label> bins;
        DynamicList<label> validBins(Pstream::nProcs());
        bin
        (
            offsets(),
            globalIds,
            order,
            bins,
            validBins
        );

        // Send local indices to individual processors as local index
        PstreamBuffers sendBufs(Pstream::commsTypes::nonBlocking, tag, comm);

        for (const auto proci : validBins)
        {
            const labelUList& es = bins[proci];

            labelList localIDs(es.size());
            forAll(es, i)
            {
                localIDs[i] = toLocal(proci, es[i]);
            }

            UOPstream os(proci, sendBufs);
            os << localIDs;
        }
        labelList recvSizes;
        sendBufs.finishedSends(recvSizes);


        PstreamBuffers returnBufs(Pstream::commsTypes::nonBlocking, tag, comm);

        forAll(recvSizes, proci)
        {
            if (recvSizes[proci])
            {
                UIPstream is(proci, sendBufs);
                labelList localIDs(is);

                // Collect entries
                List<Type> fld(localIDs.size());
                cop(fld, localIDs);

                UOPstream os(proci, returnBufs);
                os << fld;
            }
        }
        returnBufs.finishedSends();

        // Slot back
        for (const auto proci : validBins)
        {
            label start = bins.offsets()[proci];
            const SubList<label> es
            (
                order,
                bins.offsets()[proci+1]-start,  // start
                start
            );
            UIPstream is(proci, returnBufs);
            List<Type> fld(is);

            UIndirectList<Type>(allFld, es) = fld;
        }
    }
}


// ************************************************************************* //
