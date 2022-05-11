/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "Pstream.H"
#include "PstreamBuffers.H"
#include "flipOp.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T, class CombineOp, class NegateOp>
void Foam::mapDistributeBase::flipAndCombine
(
    const labelUList& map,
    const bool hasFlip,
    const UList<T>& rhs,
    const CombineOp& cop,
    const NegateOp& negOp,
    List<T>& lhs
)
{
    if (hasFlip)
    {
        forAll(map, i)
        {
            if (map[i] > 0)
            {
                label index = map[i]-1;
                cop(lhs[index], rhs[i]);
            }
            else if (map[i] < 0)
            {
                label index = -map[i]-1;
                cop(lhs[index], negOp(rhs[i]));
            }
            else
            {
                FatalErrorInFunction
                    << "At index " << i << " out of " << map.size()
                    << " have illegal index " << map[i]
                    << " for field " << rhs.size() << " with flipMap"
                    << exit(FatalError);
            }
        }
    }
    else
    {
        forAll(map, i)
        {
            cop(lhs[map[i]], rhs[i]);
        }
    }
}


template<class T, class NegateOp>
T Foam::mapDistributeBase::accessAndFlip
(
    const UList<T>& values,
    const label index,
    const bool hasFlip,
    const NegateOp& negOp
)
{
    if (hasFlip)
    {
        if (index > 0)
        {
            return values[index-1];
        }
        else if (index < 0)
        {
            return negOp(values[-index-1]);
        }
        else
        {
            FatalErrorInFunction
                << "Illegal index " << index
                << " into field of size " << values.size()
                << " with face-flipping"
                << exit(FatalError);
        }
    }

    return values[index];
}


template<class T, class NegateOp>
Foam::List<T> Foam::mapDistributeBase::accessAndFlip
(
    const UList<T>& values,
    const labelUList& indices,
    const bool hasFlip,
    const NegateOp& negOp
)
{
    const label len = indices.size();

    List<T> output(len);

    if (hasFlip)
    {
        for (label i = 0; i < len; ++i)
        {
            const label index = indices[i];

            if (index > 0)
            {
                output[i] = values[index-1];
            }
            else if (index < 0)
            {
                output[i] = negOp(values[-index-1]);
            }
            else
            {
                FatalErrorInFunction
                    << "Illegal index " << index
                    << " into field of size " << values.size()
                    << " with flipping"
                    << exit(FatalError);
            }
        }
    }
    else
    {
        // Like indirect list
        for (label i = 0; i < len; ++i)
        {
            output[i] = values[indices[i]];
        }
    }

    return output;
}


template<class T, class NegateOp>
void Foam::mapDistributeBase::distribute
(
    const Pstream::commsTypes commsType,
    const List<labelPair>& schedule,
    const label constructSize,
    const labelListList& subMap,
    const bool subHasFlip,
    const labelListList& constructMap,
    const bool constructHasFlip,
    List<T>& field,
    const NegateOp& negOp,
    const int tag,
    const label comm
)
{
    const label myRank = Pstream::myProcNo(comm);
    const label nProcs = Pstream::nProcs(comm);

    if (!Pstream::parRun())
    {
        // Do only me to me.

        List<T> subField
        (
            accessAndFlip(field, subMap[myRank], subHasFlip, negOp)
        );

        // Receive sub field from myself (subField)
        const labelList& map = constructMap[myRank];

        field.setSize(constructSize);

        flipAndCombine
        (
            map,
            constructHasFlip,
            subField,
            eqOp<T>(),
            negOp,
            field
        );

        return;
    }

    if (commsType == Pstream::commsTypes::blocking)
    {
        // Since buffered sending can reuse the field to collect the
        // received data.

        // Send sub field to neighbour
        for (const int domain : Pstream::allProcs(comm))
        {
            const labelList& map = subMap[domain];

            if (domain != myRank && map.size())
            {
                OPstream toNbr
                (
                    Pstream::commsTypes::blocking,
                    domain,
                    0,
                    tag,
                    comm
                );

                List<T> subField
                (
                    accessAndFlip(field, map, subHasFlip, negOp)
                );

                toNbr << subField;
            }
        }

        {
            // Subset myself
            List<T> subField
            (
                accessAndFlip(field, subMap[myRank], subHasFlip, negOp)
            );

            // Receive sub field from myself (subField)
            const labelList& map = constructMap[myRank];

            field.setSize(constructSize);

            flipAndCombine
            (
                map,
                constructHasFlip,
                subField,
                eqOp<T>(),
                negOp,
                field
            );
        }

        // Receive sub field from neighbour
        for (const int domain : Pstream::allProcs(comm))
        {
            const labelList& map = constructMap[domain];

            if (domain != myRank && map.size())
            {
                IPstream fromNbr
                (
                    Pstream::commsTypes::blocking,
                    domain,
                    0,
                    tag,
                    comm
                );
                List<T> subField(fromNbr);

                checkReceivedSize(domain, map.size(), subField.size());

                flipAndCombine
                (
                    map,
                    constructHasFlip,
                    subField,
                    eqOp<T>(),
                    negOp,
                    field
                );
            }
        }
    }
    else if (commsType == Pstream::commsTypes::scheduled)
    {
        // Need to make sure I don't overwrite field with received data
        // since the data might need to be sent to another processor. So
        // allocate a new field for the results.
        List<T> newField(constructSize);

        // Receive sub field from myself
        {
            List<T> subField
            (
                accessAndFlip(field, subMap[myRank], subHasFlip, negOp)
            );

            // Receive sub field from myself (subField)
            const labelList& map = constructMap[myRank];

            flipAndCombine
            (
                map,
                constructHasFlip,
                subField,
                eqOp<T>(),
                negOp,
                newField
            );
        }

        // Schedule will already have pruned 0-sized comms
        for (const labelPair& twoProcs : schedule)
        {
            // twoProcs is a swap pair of processors. The first one is the
            // one that needs to send first and then receive.

            const label sendProc = twoProcs[0];
            const label recvProc = twoProcs[1];

            if (myRank == sendProc)
            {
                // I am send first, receive next
                {
                    OPstream toNbr
                    (
                        Pstream::commsTypes::scheduled,
                        recvProc,
                        0,
                        tag,
                        comm
                    );

                    const labelList& map = subMap[recvProc];
                    List<T> subField
                    (
                        accessAndFlip(field, map, subHasFlip, negOp)
                    );

                    toNbr << subField;
                }
                {
                    IPstream fromNbr
                    (
                        Pstream::commsTypes::scheduled,
                        recvProc,
                        0,
                        tag,
                        comm
                    );
                    List<T> subField(fromNbr);

                    const labelList& map = constructMap[recvProc];

                    checkReceivedSize(recvProc, map.size(), subField.size());

                    flipAndCombine
                    (
                        map,
                        constructHasFlip,
                        subField,
                        eqOp<T>(),
                        negOp,
                        newField
                    );
                }
            }
            else
            {
                // I am receive first, send next
                {
                    IPstream fromNbr
                    (
                        Pstream::commsTypes::scheduled,
                        sendProc,
                        0,
                        tag,
                        comm
                    );
                    List<T> subField(fromNbr);

                    const labelList& map = constructMap[sendProc];

                    checkReceivedSize(sendProc, map.size(), subField.size());

                    flipAndCombine
                    (
                        map,
                        constructHasFlip,
                        subField,
                        eqOp<T>(),
                        negOp,
                        newField
                    );
                }
                {
                    OPstream toNbr
                    (
                        Pstream::commsTypes::scheduled,
                        sendProc,
                        0,
                        tag,
                        comm
                    );

                    const labelList& map = subMap[sendProc];
                    List<T> subField
                    (
                        accessAndFlip(field, map, subHasFlip, negOp)
                    );

                    toNbr << subField;
                }
            }
        }
        field.transfer(newField);
    }
    else if (commsType == Pstream::commsTypes::nonBlocking)
    {
        const label nOutstanding = Pstream::nRequests();

        if (!is_contiguous<T>::value)
        {
            PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking, tag, comm);

            // Stream data into buffer
            for (const int domain : Pstream::allProcs(comm))
            {
                const labelList& map = subMap[domain];

                if (domain != myRank && map.size())
                {
                    // Put data into send buffer
                    UOPstream toDomain(domain, pBufs);

                    List<T> subField
                    (
                        accessAndFlip(field, map, subHasFlip, negOp)
                    );

                    toDomain << subField;
                }
            }

            // Start receiving. Do not block.
            pBufs.finishedSends(false);

            {
                // Set up 'send' to myself
                List<T> mySubField
                (
                    accessAndFlip(field, subMap[myRank], subHasFlip, negOp)
                );

                // Combine bits. Note that can reuse field storage
                field.setSize(constructSize);

                // Receive sub field from myself
                const labelList& map = constructMap[myRank];

                flipAndCombine
                (
                    map,
                    constructHasFlip,
                    mySubField,
                    eqOp<T>(),
                    negOp,
                    field
                );
            }

            // Block ourselves, waiting only for the current comms
            Pstream::waitRequests(nOutstanding);

            // Consume
            for (const int domain : Pstream::allProcs(comm))
            {
                const labelList& map = constructMap[domain];

                if (domain != myRank && map.size())
                {
                    UIPstream str(domain, pBufs);
                    List<T> recvField(str);

                    checkReceivedSize(domain, map.size(), recvField.size());

                    flipAndCombine
                    (
                        map,
                        constructHasFlip,
                        recvField,
                        eqOp<T>(),
                        negOp,
                        field
                    );
                }
            }
        }
        else
        {
            // Set up sends to neighbours

            List<List<T>> sendFields(nProcs);

            for (const int domain : Pstream::allProcs(comm))
            {
                const labelList& map = subMap[domain];

                if (domain != myRank && map.size())
                {
                    sendFields[domain] =
                        accessAndFlip(field, map, subHasFlip, negOp);

                    OPstream::write
                    (
                        Pstream::commsTypes::nonBlocking,
                        domain,
                        sendFields[domain].cdata_bytes(),
                        sendFields[domain].size_bytes(),
                        tag,
                        comm
                    );
                }
            }

            // Set up receives from neighbours

            List<List<T>> recvFields(nProcs);

            for (const int domain : Pstream::allProcs(comm))
            {
                const labelList& map = constructMap[domain];

                if (domain != myRank && map.size())
                {
                    recvFields[domain].resize(map.size());

                    IPstream::read
                    (
                        Pstream::commsTypes::nonBlocking,
                        domain,
                        recvFields[domain].data_bytes(),
                        recvFields[domain].size_bytes(),
                        tag,
                        comm
                    );
                }
            }


            // Set up 'send' to myself
            {
                sendFields[myRank] =
                    accessAndFlip(field, subMap[myRank], subHasFlip, negOp);
            }


            // Combine bits. Note that can reuse field storage

            field.setSize(constructSize);


            // Receive sub field from myself (sendFields[myRank])
            {
                const labelList& map = constructMap[myRank];
                const List<T>& subField = sendFields[myRank];

                flipAndCombine
                (
                    map,
                    constructHasFlip,
                    subField,
                    eqOp<T>(),
                    negOp,
                    field
                );
            }


            // Wait for all to finish

            Pstream::waitRequests(nOutstanding);


            // Collect neighbour fields

            for (const int domain : Pstream::allProcs(comm))
            {
                const labelList& map = constructMap[domain];

                if (domain != myRank && map.size())
                {
                    const List<T>& subField = recvFields[domain];

                    checkReceivedSize(domain, map.size(), subField.size());

                    flipAndCombine
                    (
                        map,
                        constructHasFlip,
                        subField,
                        eqOp<T>(),
                        negOp,
                        field
                    );
                }
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unknown communication schedule " << int(commsType)
            << abort(FatalError);
    }
}


template<class T, class CombineOp, class NegateOp>
void Foam::mapDistributeBase::distribute
(
    const Pstream::commsTypes commsType,
    const List<labelPair>& schedule,
    const label constructSize,
    const labelListList& subMap,
    const bool subHasFlip,
    const labelListList& constructMap,
    const bool constructHasFlip,
    List<T>& field,
    const T& nullValue,
    const CombineOp& cop,
    const NegateOp& negOp,
    const int tag,
    const label comm
)
{
    const label myRank = Pstream::myProcNo(comm);
    const label nProcs = Pstream::nProcs(comm);

    if (!Pstream::parRun())
    {
        // Do only me to me.

        List<T> subField
        (
            accessAndFlip(field, subMap[myRank], subHasFlip, negOp)
        );

        // Receive sub field from myself (subField)
        const labelList& map = constructMap[myRank];

        field.resize_nocopy(constructSize);
        field = nullValue;

        flipAndCombine(map, constructHasFlip, subField, cop, negOp, field);

        return;
    }

    if (commsType == Pstream::commsTypes::blocking)
    {
        // Since buffered sending can reuse the field to collect the
        // received data.

        // Send sub field to neighbour
        for (const int domain : Pstream::allProcs(comm))
        {
            const labelList& map = subMap[domain];

            if (domain != myRank && map.size())
            {
                OPstream toNbr
                (
                    Pstream::commsTypes::blocking,
                    domain,
                    0,
                    tag,
                    comm
                );
                List<T> subField
                (
                    accessAndFlip(field, map, subHasFlip, negOp)
                );

                toNbr << subField;
            }
        }

        {
            // Subset myself
            List<T> subField
            (
                accessAndFlip(field, subMap[myRank], subHasFlip, negOp)
            );

            // Receive sub field from myself (subField)
            const labelList& map = constructMap[myRank];

            field.resize_nocopy(constructSize);
            field = nullValue;

            flipAndCombine
            (
                map,
                constructHasFlip,
                subField,
                cop,
                negOp,
                field
            );
        }

        // Receive sub field from neighbour
        for (const int domain : Pstream::allProcs(comm))
        {
            const labelList& map = constructMap[domain];

            if (domain != myRank && map.size())
            {
                IPstream fromNbr
                (
                    Pstream::commsTypes::blocking,
                    domain,
                    0,
                    tag,
                    comm
                );
                List<T> subField(fromNbr);

                checkReceivedSize(domain, map.size(), subField.size());

                flipAndCombine
                (
                    map,
                    constructHasFlip,
                    subField,
                    cop,
                    negOp,
                    field
                );
            }
        }
    }
    else if (commsType == Pstream::commsTypes::scheduled)
    {
        // Need to make sure I don't overwrite field with received data
        // since the data might need to be sent to another processor. So
        // allocate a new field for the results.
        List<T> newField(constructSize, nullValue);

        {
            // Subset myself
            List<T> subField
            (
                accessAndFlip(field, subMap[myRank], subHasFlip, negOp)
            );

            // Receive sub field from myself (subField)
            const labelList& map = constructMap[myRank];

            flipAndCombine
            (
                map,
                constructHasFlip,
                subField,
                cop,
                negOp,
                newField
            );
        }


        // Schedule will already have pruned 0-sized comms
        for (const labelPair& twoProcs : schedule)
        {
            // twoProcs is a swap pair of processors. The first one is the
            // one that needs to send first and then receive.

            const label sendProc = twoProcs[0];
            const label recvProc = twoProcs[1];

            if (myRank == sendProc)
            {
                // I am send first, receive next
                {
                    OPstream toNbr
                    (
                        Pstream::commsTypes::scheduled,
                        recvProc,
                        0,
                        tag,
                        comm
                    );

                    const labelList& map = subMap[recvProc];

                    List<T> subField
                    (
                        accessAndFlip(field, map, subHasFlip, negOp)
                    );

                    toNbr << subField;
                }
                {
                    IPstream fromNbr
                    (
                        Pstream::commsTypes::scheduled,
                        recvProc,
                        0,
                        tag,
                        comm
                    );
                    List<T> subField(fromNbr);
                    const labelList& map = constructMap[recvProc];

                    checkReceivedSize(recvProc, map.size(), subField.size());

                    flipAndCombine
                    (
                        map,
                        constructHasFlip,
                        subField,
                        cop,
                        negOp,
                        newField
                    );
                }
            }
            else
            {
                // I am receive first, send next
                {
                    IPstream fromNbr
                    (
                        Pstream::commsTypes::scheduled,
                        sendProc,
                        0,
                        tag,
                        comm
                    );
                    List<T> subField(fromNbr);
                    const labelList& map = constructMap[sendProc];

                    checkReceivedSize(sendProc, map.size(), subField.size());

                    flipAndCombine
                    (
                        map,
                        constructHasFlip,
                        subField,
                        cop,
                        negOp,
                        newField
                    );
                }
                {
                    OPstream toNbr
                    (
                        Pstream::commsTypes::scheduled,
                        sendProc,
                        0,
                        tag,
                        comm
                    );

                    const labelList& map = subMap[sendProc];

                    List<T> subField
                    (
                        accessAndFlip(field, map, subHasFlip, negOp)
                    );

                    toNbr << subField;
                }
            }
        }
        field.transfer(newField);
    }
    else if (commsType == Pstream::commsTypes::nonBlocking)
    {
        const label nOutstanding = Pstream::nRequests();

        if (!is_contiguous<T>::value)
        {
            PstreamBuffers pBufs(Pstream::commsTypes::nonBlocking, tag, comm);

            // Stream data into buffer
            for (const int domain : Pstream::allProcs(comm))
            {
                const labelList& map = subMap[domain];

                if (domain != myRank && map.size())
                {
                    // Put data into send buffer
                    UOPstream toDomain(domain, pBufs);

                    List<T> subField
                    (
                        accessAndFlip(field, map, subHasFlip, negOp)
                    );

                    toDomain << subField;
                }
            }

            // Start receiving. Do not block.
            pBufs.finishedSends(false);

            {
                // Set up 'send' to myself
                List<T> mySubField
                (
                    accessAndFlip(field, subMap[myRank], subHasFlip, negOp)
                );

                // Combine bits. Note that can reuse field storage
                field.resize_nocopy(constructSize);
                field = nullValue;

                // Receive sub field from myself
                const labelList& map = constructMap[myRank];

                flipAndCombine
                (
                    map,
                    constructHasFlip,
                    mySubField,
                    cop,
                    negOp,
                    field
                );
            }

            // Block ourselves, waiting only for the current comms
            Pstream::waitRequests(nOutstanding);

            // Consume
            for (const int domain : Pstream::allProcs(comm))
            {
                const labelList& map = constructMap[domain];

                if (domain != myRank && map.size())
                {
                    UIPstream str(domain, pBufs);
                    List<T> recvField(str);

                    checkReceivedSize(domain, map.size(), recvField.size());

                    flipAndCombine
                    (
                        map,
                        constructHasFlip,
                        recvField,
                        cop,
                        negOp,
                        field
                    );
                }
            }
        }
        else
        {
            // Set up sends to neighbours

            List<List<T>> sendFields(nProcs);

            for (const int domain : Pstream::allProcs(comm))
            {
                const labelList& map = subMap[domain];

                if (domain != myRank && map.size())
                {
                    sendFields[domain] =
                        accessAndFlip(field, map, subHasFlip, negOp);

                    OPstream::write
                    (
                        Pstream::commsTypes::nonBlocking,
                        domain,
                        sendFields[domain].cdata_bytes(),
                        sendFields[domain].size_bytes(),
                        tag,
                        comm
                    );
                }
            }

            // Set up receives from neighbours

            List<List<T>> recvFields(nProcs);

            for (const int domain : Pstream::allProcs(comm))
            {
                const labelList& map = constructMap[domain];

                if (domain != myRank && map.size())
                {
                    recvFields[domain].setSize(map.size());
                    UIPstream::read
                    (
                        Pstream::commsTypes::nonBlocking,
                        domain,
                        recvFields[domain].data_bytes(),
                        recvFields[domain].size_bytes(),
                        tag,
                        comm
                    );
                }
            }

            // Set up 'send' to myself

            {
                sendFields[myRank] =
                    accessAndFlip(field, subMap[myRank], subHasFlip, negOp);
            }


            // Combine bits. Note that can reuse field storage

            field.resize_nocopy(constructSize);
            field = nullValue;

            // Receive sub field from myself (subField)
            {
                const labelList& map = constructMap[myRank];
                const List<T>& subField = sendFields[myRank];

                flipAndCombine
                (
                    map,
                    constructHasFlip,
                    subField,
                    cop,
                    negOp,
                    field
                );
            }


            // Wait for all to finish

            Pstream::waitRequests(nOutstanding);


            // Collect neighbour fields

            for (const int domain : Pstream::allProcs(comm))
            {
                const labelList& map = constructMap[domain];

                if (domain != myRank && map.size())
                {
                    const List<T>& subField = recvFields[domain];

                    checkReceivedSize(domain, map.size(), subField.size());

                    flipAndCombine
                    (
                        map,
                        constructHasFlip,
                        subField,
                        cop,
                        negOp,
                        field
                    );
                }
            }
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unknown communication schedule " << int(commsType)
            << abort(FatalError);
    }
}


template<class T>
void Foam::mapDistributeBase::send(PstreamBuffers& pBufs, const List<T>& field)
const
{
    // Stream data into buffer
    for (const int domain : Pstream::allProcs(comm_))
    {
        const labelList& map = subMap_[domain];

        if (map.size())
        {
            // Put data into send buffer
            UOPstream toDomain(domain, pBufs);

            List<T> subField
            (
                accessAndFlip(field, map, subHasFlip_, flipOp())
            );

            toDomain << subField;
        }
    }

    // Start sending and receiving but do not block.
    pBufs.finishedSends(false);
}


template<class T>
void Foam::mapDistributeBase::receive(PstreamBuffers& pBufs, List<T>& field)
const
{
    // Consume
    field.resize_nocopy(constructSize_);

    for (const int domain : Pstream::allProcs(comm_))
    {
        const labelList& map = constructMap_[domain];

        if (map.size())
        {
            UIPstream str(domain, pBufs);
            List<T> recvField(str);

            if (recvField.size() != map.size())
            {
                FatalErrorInFunction
                    << "Expected from processor " << domain
                    << " " << map.size() << " but received "
                    << recvField.size() << " elements."
                    << abort(FatalError);
            }

            flipAndCombine
            (
                map,
                constructHasFlip_,
                recvField,
                eqOp<T>(),
                flipOp(),
                field
            );
        }
    }
}


template<class T, class NegateOp>
void Foam::mapDistributeBase::distribute
(
    const Pstream::commsTypes commsType,
    List<T>& values,
    const NegateOp& negOp,
    const int tag
) const
{
    distribute
    (
        commsType,
        whichSchedule(commsType),
        constructSize_,
        subMap_,
        subHasFlip_,
        constructMap_,
        constructHasFlip_,
        values,
        negOp,
        tag,
        comm_
    );
}


template<class T, class NegateOp>
void Foam::mapDistributeBase::distribute
(
    const Pstream::commsTypes commsType,
    const T& nullValue,
    List<T>& values,
    const NegateOp& negOp,
    const int tag
) const
{
    distribute
    (
        commsType,
        whichSchedule(commsType),
        constructSize_,
        subMap_,
        subHasFlip_,
        constructMap_,
        constructHasFlip_,
        values,
        nullValue,
        eqOp<T>(),
        negOp,
        tag,
        comm_
    );
}


template<class T, class NegateOp>
void Foam::mapDistributeBase::distribute
(
    List<T>& values,
    const NegateOp& negOp,
    const int tag
) const
{
    distribute
    (
        UPstream::defaultCommsType, values, negOp, tag
    );
}


template<class T>
void Foam::mapDistributeBase::distribute
(
    List<T>& values,
    const int tag
) const
{
    distribute(values, flipOp(), tag);
}


template<class T>
void Foam::mapDistributeBase::distribute
(
    DynamicList<T>& values,
    const int tag
) const
{
    values.shrink();

    List<T>& list = static_cast<List<T>&>(values);

    distribute(list, tag);

    values.setCapacity(list.size());
}


template<class T>
void Foam::mapDistributeBase::reverseDistribute
(
    const Pstream::commsTypes commsType,
    const label constructSize,
    List<T>& values,
    const int tag
) const
{
    reverseDistribute<T, flipOp>
    (
        commsType,
        constructSize,
        values,
        flipOp(),
        tag
    );
}


template<class T, class NegateOp>
void Foam::mapDistributeBase::reverseDistribute
(
    const Pstream::commsTypes commsType,
    const label constructSize,
    List<T>& values,
    const NegateOp& negOp,
    const int tag
) const
{
    distribute
    (
        commsType,
        whichSchedule(commsType),
        constructSize,
        constructMap_,
        constructHasFlip_,
        subMap_,
        subHasFlip_,
        values,
        negOp,
        tag,
        comm_
    );
}


template<class T>
void Foam::mapDistributeBase::reverseDistribute
(
    const Pstream::commsTypes commsType,
    const label constructSize,
    const T& nullValue,
    List<T>& values,
    const int tag
) const
{
    distribute
    (
        commsType,
        whichSchedule(commsType),
        constructSize,
        constructMap_,
        constructHasFlip_,
        subMap_,
        subHasFlip_,
        values,

        nullValue,
        eqOp<T>(),
        flipOp(),

        tag,
        comm_
    );
}


template<class T>
void Foam::mapDistributeBase::reverseDistribute
(
    const label constructSize,
    List<T>& values,
    const int tag
) const
{
    reverseDistribute
    (
        UPstream::defaultCommsType,
        constructSize,
        values,
        tag
    );
}


template<class T>
void Foam::mapDistributeBase::reverseDistribute
(
    const label constructSize,
    const T& nullValue,
    List<T>& values,
    const int tag
) const
{
    reverseDistribute
    (
        UPstream::defaultCommsType,
        constructSize,
        nullValue,
        values,
        tag
    );
}


// ************************************************************************* //
