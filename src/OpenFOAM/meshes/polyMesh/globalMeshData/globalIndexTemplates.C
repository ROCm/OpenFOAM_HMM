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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Container, class Type>
void Foam::globalIndex::gather
(
    const labelUList& off,
    const label comm,
    const Container& procIDs,
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
        allFld.resize(off.last());

        // Assign my local data
        SubList<Type>(allFld, fld.size(), 0) = fld;

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
                        reinterpret_cast<char*>(procSlot.data()),
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
                    reinterpret_cast<char*>(procSlot.data()),
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
                    reinterpret_cast<const char*>(fld.cdata()),
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
                reinterpret_cast<const char*>(fld.cdata()),
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
        allFld.resize(off.last());

        // Assign my local data
        SubList<Type>(allFld, fld.size(), 0) = fld;

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
    const Pstream::commsTypes commsType
) const
{
    gather
    (
        UPstream::worldComm,
        UPstream::procID(UPstream::worldComm),
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
    const Pstream::commsTypes commsType
) const
{
    gather
    (
        offsets_,
        UPstream::worldComm,
        UPstream::procID(UPstream::worldComm),
        fld,
        allFld,
        tag,
        commsType
    );
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


template<class Container, class Type>
void Foam::globalIndex::gather
(
    const labelUList& off,
    const label comm,
    const Container& procIDs,
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
    const Pstream::commsTypes commsType
) const
{
    List<Type> allFld;

    gather
    (
        UPstream::worldComm,
        UPstream::procID(UPstream::worldComm),
        fld,
        allFld,
        tag,
        commsType
    );

    if (Pstream::master(UPstream::worldComm))
    {
        fld.transfer(allFld);
    }
    else
    {
        fld.clear();
    }
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


template<class Container, class Type>
void Foam::globalIndex::scatter
(
    const labelUList& off,
    const label comm,
    const Container& procIDs,
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
        fld.deepCopy(SubList<Type>(allFld, off[1]-off[0]));

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
                        reinterpret_cast<const char*>(procSlot.cdata()),
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
                    reinterpret_cast<const char*>(procSlot.cdata()),
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
                    reinterpret_cast<char*>(fld.data()),
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
                reinterpret_cast<char*>(fld.data()),
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
    const Pstream::commsTypes commsType
) const
{
    scatter
    (
        offsets_,
        UPstream::worldComm,
        UPstream::procID(UPstream::worldComm),
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
    allFld.resize(globalIds.size());
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
