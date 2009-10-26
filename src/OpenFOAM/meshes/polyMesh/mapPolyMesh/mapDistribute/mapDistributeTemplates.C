/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "Pstream.H"
#include "PstreamBuffers.H"
#include "PstreamCombineReduceOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Distribute list.
template<class T>
void Foam::mapDistribute::distribute
(
    const Pstream::commsTypes commsType,
    const List<labelPair>& schedule,
    const label constructSize,
    const labelListList& subMap,
    const labelListList& constructMap,
    List<T>& field
)
{
    if (commsType == Pstream::blocking)
    {
        // Since buffered sending can reuse the field to collect the
        // received data.

        // Send sub field to neighbour
        for (label domain = 0; domain < Pstream::nProcs(); domain++)
        {
            const labelList& map = subMap[domain];

            if (domain != Pstream::myProcNo() && map.size())
            {
                OPstream toNbr(Pstream::blocking, domain);
                toNbr << UIndirectList<T>(field, map);
            }
        }

        // Subset myself
        const labelList& mySubMap = subMap[Pstream::myProcNo()];

        List<T> subField(mySubMap.size());
        forAll(mySubMap, i)
        {
            subField[i] = field[mySubMap[i]];
        }

        // Receive sub field from myself (subField)
        const labelList& map = constructMap[Pstream::myProcNo()];

        field.setSize(constructSize);

        forAll(map, i)
        {
            field[map[i]] = subField[i];
        }

        // Receive sub field from neighbour
        for (label domain = 0; domain < Pstream::nProcs(); domain++)
        {
            const labelList& map = constructMap[domain];

            if (domain != Pstream::myProcNo() && map.size())
            {
                IPstream fromNbr(Pstream::blocking, domain);
                List<T> subField(fromNbr);

                if (subField.size() != map.size())
                {
                    FatalErrorIn
                    (
                        "template<class T>\n"
                        "void mapDistribute::distribute\n"
                        "(\n"
                        "    const Pstream::commsTypes commsType,\n"
                        "    const List<labelPair>& schedule,\n"
                        "    const label constructSize,\n"
                        "    const labelListList& subMap,\n"
                        "    const labelListList& constructMap,\n"
                        "    List<T>& field\n"
                        ")\n"
                    )   << "Expected from processor " << domain
                        << " " << map.size() << " but received "
                        << subField.size() << " elements."
                        << abort(FatalError);
                }

                forAll(map, i)
                {
                    field[map[i]] = subField[i];
                }
            }
        }
    }
    else if (commsType == Pstream::scheduled)
    {
        // Need to make sure I don't overwrite field with received data
        // since the data might need to be sent to another processor. So
        // allocate a new field for the results.
        List<T> newField(constructSize);

        // Subset myself
        UIndirectList<T> subField(field, subMap[Pstream::myProcNo()]);

        // Receive sub field from myself (subField)
        const labelList& map = constructMap[Pstream::myProcNo()];

        forAll(map, i)
        {
            newField[map[i]] = subField[i];
        }

        // Schedule will already have pruned 0-sized comms
        forAll(schedule, i)
        {
            const labelPair& twoProcs = schedule[i];
            label sendProc = twoProcs[0];
            label recvProc = twoProcs[1];

            if (Pstream::myProcNo() == sendProc)
            {
                // I am sender. Send to recvProc.
                OPstream toNbr(Pstream::scheduled, recvProc);
                toNbr << UIndirectList<T>(field, subMap[recvProc]);
            }
            else
            {
                // I am receiver. Receive from sendProc.
                IPstream fromNbr(Pstream::scheduled, sendProc);
                List<T> subField(fromNbr);

                const labelList& map = constructMap[sendProc];

                if (subField.size() != map.size())
                {
                    FatalErrorIn
                    (
                        "template<class T>\n"
                        "void mapDistribute::distribute\n"
                        "(\n"
                        "    const Pstream::commsTypes commsType,\n"
                        "    const List<labelPair>& schedule,\n"
                        "    const label constructSize,\n"
                        "    const labelListList& subMap,\n"
                        "    const labelListList& constructMap,\n"
                        "    List<T>& field\n"
                        ")\n"
                    )   << "Expected from processor " << sendProc
                        << " " << map.size() << " but received "
                        << subField.size() << " elements."
                        << abort(FatalError);
                }

                forAll(map, i)
                {
                    newField[map[i]] = subField[i];
                }
            }
        }
        field.transfer(newField);
    }
    else if (commsType == Pstream::nonBlocking)
    {
        if (!contiguous<T>())
        {
            PstreamBuffers pBuffs(Pstream::nonBlocking);

            // Stream data into buffer
            for (label domain = 0; domain < Pstream::nProcs(); domain++)
            {
                const labelList& map = subMap[domain];

                if (domain != Pstream::myProcNo() && map.size())
                {
                    // Put data into send buffer
                    UOPstream toDomain(domain, pBuffs);
                    toDomain << UIndirectList<T>(field, map);
                }
            }

            // Start receiving
            pBuffs.finishedSends();

            {
                // Set up 'send' to myself
                const labelList& mySubMap = subMap[Pstream::myProcNo()];
                List<T> mySubField(mySubMap.size());
                forAll(mySubMap, i)
                {
                    mySubField[i] = field[mySubMap[i]];
                }
                // Combine bits. Note that can reuse field storage
                field.setSize(constructSize);
                // Receive sub field from myself
                {
                    const labelList& map = constructMap[Pstream::myProcNo()];

                    forAll(map, i)
                    {
                        field[map[i]] = mySubField[i];
                    }
                }
            }

            // Consume
            for (label domain = 0; domain < Pstream::nProcs(); domain++)
            {
                const labelList& map = constructMap[domain];

                if (domain != Pstream::myProcNo() && map.size())
                {
                    UIPstream str(domain, pBuffs);
                    List<T> recvField(str);

                    if (recvField.size() != map.size())
                    {
                        FatalErrorIn
                        (
                            "template<class T>\n"
                            "void mapDistribute::distribute\n"
                            "(\n"
                            "    const Pstream::commsTypes commsType,\n"
                            "    const List<labelPair>& schedule,\n"
                            "    const label constructSize,\n"
                            "    const labelListList& subMap,\n"
                            "    const labelListList& constructMap,\n"
                            "    List<T>& field\n"
                            ")\n"
                        )   << "Expected from processor " << domain
                            << " " << map.size() << " but received "
                            << recvField.size() << " elements."
                            << abort(FatalError);
                    }

                    forAll(map, i)
                    {
                        field[map[i]] = recvField[i];
                    }
                }
            }
        }
        else
        {
            // Set up sends to neighbours

            List<List<T > > sendFields(Pstream::nProcs());

            for (label domain = 0; domain < Pstream::nProcs(); domain++)
            {
                const labelList& map = subMap[domain];

                if (domain != Pstream::myProcNo() && map.size())
                {
                    List<T>& subField = sendFields[domain];
                    subField.setSize(map.size());
                    forAll(map, i)
                    {
                        subField[i] = field[map[i]];
                    }

                    OPstream::write
                    (
                        Pstream::nonBlocking,
                        domain,
                        reinterpret_cast<const char*>(subField.begin()),
                        subField.byteSize()
                    );
                }
            }

            // Set up receives from neighbours

            List<List<T > > recvFields(Pstream::nProcs());

            for (label domain = 0; domain < Pstream::nProcs(); domain++)
            {
                const labelList& map = constructMap[domain];

                if (domain != Pstream::myProcNo() && map.size())
                {
                    recvFields[domain].setSize(map.size());
                    IPstream::read
                    (
                        Pstream::nonBlocking,
                        domain,
                        reinterpret_cast<char*>(recvFields[domain].begin()),
                        recvFields[domain].byteSize()
                    );
                }
            }


            // Set up 'send' to myself

            {
                const labelList& map = subMap[Pstream::myProcNo()];

                List<T>& subField = sendFields[Pstream::myProcNo()];
                subField.setSize(map.size());
                forAll(map, i)
                {
                    subField[i] = field[map[i]];
                }
            }


            // Combine bits. Note that can reuse field storage

            field.setSize(constructSize);


            // Receive sub field from myself (sendFields[Pstream::myProcNo()])
            {
                const labelList& map = constructMap[Pstream::myProcNo()];
                const List<T>& subField = sendFields[Pstream::myProcNo()];

                forAll(map, i)
                {
                    field[map[i]] = subField[i];
                }
            }


            // Wait for all to finish

            Pstream::waitRequests();

            // Collect neighbour fields

            for (label domain = 0; domain < Pstream::nProcs(); domain++)
            {
                const labelList& map = constructMap[domain];

                if (domain != Pstream::myProcNo() && map.size())
                {
                    if (recvFields[domain].size() != map.size())
                    {
                        FatalErrorIn
                        (
                            "template<class T>\n"
                            "void mapDistribute::distribute\n"
                            "(\n"
                            "    const Pstream::commsTypes commsType,\n"
                            "    const List<labelPair>& schedule,\n"
                            "    const label constructSize,\n"
                            "    const labelListList& subMap,\n"
                            "    const labelListList& constructMap,\n"
                            "    List<T>& field\n"
                            ")\n"
                        )   << "Expected from processor " << domain
                            << " " << map.size() << " but received "
                            << recvFields[domain].size() << " elements."
                            << abort(FatalError);
                    }

                    forAll(map, i)
                    {
                        field[map[i]] = recvFields[domain][i];
                    }
                }
            }
        }
    }
    else
    {
        FatalErrorIn("mapDistribute::distribute(..)")
            << "Unknown communication schedule " << commsType
            << abort(FatalError);
    }
}


// Distribute list.
template<class T, class CombineOp>
void Foam::mapDistribute::distribute
(
    const Pstream::commsTypes commsType,
    const List<labelPair>& schedule,
    const label constructSize,
    const labelListList& subMap,
    const labelListList& constructMap,
    List<T>& field,
    const CombineOp& cop,
    const T& nullValue
)
{
    if (commsType == Pstream::blocking)
    {
        // Since buffered sending can reuse the field to collect the
        // received data.

        // Send sub field to neighbour
        for (label domain = 0; domain < Pstream::nProcs(); domain++)
        {
            const labelList& map = subMap[domain];

            if (domain != Pstream::myProcNo() && map.size())
            {
                OPstream toNbr(Pstream::blocking, domain);
                toNbr << UIndirectList<T>(field, map);
            }
        }

        // Subset myself
        const labelList& mySubMap = subMap[Pstream::myProcNo()];

        List<T> subField(mySubMap.size());
        forAll(mySubMap, i)
        {
            subField[i] = field[mySubMap[i]];
        }

        // Receive sub field from myself (subField)
        const labelList& map = constructMap[Pstream::myProcNo()];

        field.setSize(constructSize);
        field = nullValue;

        forAll(map, i)
        {
            cop(field[map[i]], subField[i]);
        }

        // Receive sub field from neighbour
        for (label domain = 0; domain < Pstream::nProcs(); domain++)
        {
            const labelList& map = constructMap[domain];

            if (domain != Pstream::myProcNo() && map.size())
            {
                IPstream fromNbr(Pstream::blocking, domain);
                List<T> subField(fromNbr);

                if (subField.size() != map.size())
                {
                    FatalErrorIn
                    (
                        "template<class T>\n"
                        "void mapDistribute::distribute\n"
                        "(\n"
                        "    const Pstream::commsTypes commsType,\n"
                        "    const List<labelPair>& schedule,\n"
                        "    const label constructSize,\n"
                        "    const labelListList& subMap,\n"
                        "    const labelListList& constructMap,\n"
                        "    List<T>& field\n"
                        ")\n"
                    )   << "Expected from processor " << domain
                        << " " << map.size() << " but received "
                        << subField.size() << " elements."
                        << abort(FatalError);
                }

                forAll(map, i)
                {
                    cop(field[map[i]], subField[i]);
                }
            }
        }
    }
    else if (commsType == Pstream::scheduled)
    {
        // Need to make sure I don't overwrite field with received data
        // since the data might need to be sent to another processor. So
        // allocate a new field for the results.
        List<T> newField(constructSize, nullValue);

        // Subset myself
        UIndirectList<T> subField(field, subMap[Pstream::myProcNo()]);

        // Receive sub field from myself (subField)
        const labelList& map = constructMap[Pstream::myProcNo()];

        forAll(map, i)
        {
            cop(newField[map[i]], subField[i]);
        }

        // Schedule will already have pruned 0-sized comms
        forAll(schedule, i)
        {
            const labelPair& twoProcs = schedule[i];
            label sendProc = twoProcs[0];
            label recvProc = twoProcs[1];

            if (Pstream::myProcNo() == sendProc)
            {
                // I am sender. Send to recvProc.
                OPstream toNbr(Pstream::scheduled, recvProc);
                toNbr << UIndirectList<T>(field, subMap[recvProc]);
            }
            else
            {
                // I am receiver. Receive from sendProc.
                IPstream fromNbr(Pstream::scheduled, sendProc);
                List<T> subField(fromNbr);

                const labelList& map = constructMap[sendProc];

                if (subField.size() != map.size())
                {
                    FatalErrorIn
                    (
                        "template<class T>\n"
                        "void mapDistribute::distribute\n"
                        "(\n"
                        "    const Pstream::commsTypes commsType,\n"
                        "    const List<labelPair>& schedule,\n"
                        "    const label constructSize,\n"
                        "    const labelListList& subMap,\n"
                        "    const labelListList& constructMap,\n"
                        "    List<T>& field\n"
                        ")\n"
                    )   << "Expected from processor " << sendProc
                        << " " << map.size() << " but received "
                        << subField.size() << " elements."
                        << abort(FatalError);
                }

                forAll(map, i)
                {
                    cop(newField[map[i]], subField[i]);
                }
            }
        }
        field.transfer(newField);
    }
    else if (commsType == Pstream::nonBlocking)
    {
        if (!contiguous<T>())
        {
//XXXXXX
            PstreamBuffers pBuffs(Pstream::nonBlocking);

            // Stream data into buffer
            for (label domain = 0; domain < Pstream::nProcs(); domain++)
            {
                const labelList& map = subMap[domain];

                if (domain != Pstream::myProcNo() && map.size())
                {
                    // Put data into send buffer
                    UOPstream toDomain(domain, pBuffs);
                    toDomain << UIndirectList<T>(field, map);
                }
            }

            // Start receiving
            pBuffs.finishedSends();

            {
                // Set up 'send' to myself
                List<T> mySubField(field, subMap[Pstream::myProcNo()]);
                // Combine bits. Note that can reuse field storage
                field.setSize(constructSize);
                field = nullValue;
                // Receive sub field from myself
                {
                    const labelList& map = constructMap[Pstream::myProcNo()];

                    forAll(map, i)
                    {
                        cop(field[map[i]], mySubField[i]);
                    }
                }
            }


            // Wait till all finished
            UPstream::waitRequests();

            // Consume
            for (label domain = 0; domain < Pstream::nProcs(); domain++)
            {
                const labelList& map = constructMap[domain];

                if (domain != Pstream::myProcNo() && map.size())
                {
                    UIPstream str(domain, pBuffs);
                    List<T> recvField(str);

                    if (recvField.size() != map.size())
                    {
                        FatalErrorIn
                        (
                            "template<class T>\n"
                            "void mapDistribute::distribute\n"
                            "(\n"
                            "    const Pstream::commsTypes commsType,\n"
                            "    const List<labelPair>& schedule,\n"
                            "    const label constructSize,\n"
                            "    const labelListList& subMap,\n"
                            "    const labelListList& constructMap,\n"
                            "    List<T>& field\n"
                            ")\n"
                        )   << "Expected from processor " << domain
                            << " " << map.size() << " but received "
                            << recvField.size() << " elements."
                            << abort(FatalError);
                    }

                    forAll(map, i)
                    {
                        cop(field[map[i]], recvField[i]);
                    }
                }
            }
        }
        else
        {
            // Set up sends to neighbours

            List<List<T > > sendFields(Pstream::nProcs());

            for (label domain = 0; domain < Pstream::nProcs(); domain++)
            {
                const labelList& map = subMap[domain];

                if (domain != Pstream::myProcNo() && map.size())
                {
                    List<T>& subField = sendFields[domain];
                    subField.setSize(map.size());
                    forAll(map, i)
                    {
                        subField[i] = field[map[i]];
                    }

                    OPstream::write
                    (
                        Pstream::nonBlocking,
                        domain,
                        reinterpret_cast<const char*>(subField.begin()),
                        subField.size()*sizeof(T)
                    );
                }
            }

            // Set up receives from neighbours

            List<List<T > > recvFields(Pstream::nProcs());

            for (label domain = 0; domain < Pstream::nProcs(); domain++)
            {
                const labelList& map = constructMap[domain];

                if (domain != Pstream::myProcNo() && map.size())
                {
                    recvFields[domain].setSize(map.size());
                    UIPstream::read
                    (
                        Pstream::nonBlocking,
                        domain,
                        reinterpret_cast<char*>(recvFields[domain].begin()),
                        recvFields[domain].size()*sizeof(T)
                    );
                }
            }

            // Set up 'send' to myself

            {
                const labelList& map = subMap[Pstream::myProcNo()];

                List<T>& subField = sendFields[Pstream::myProcNo()];
                subField.setSize(map.size());
                forAll(map, i)
                {
                    subField[i] = field[map[i]];
                }
            }


            // Combine bits. Note that can reuse field storage

            field.setSize(constructSize);
            field = nullValue;

            // Receive sub field from myself (subField)
            {
                const labelList& map = constructMap[Pstream::myProcNo()];
                const List<T>& subField = sendFields[Pstream::myProcNo()];

                forAll(map, i)
                {
                    cop(field[map[i]], subField[i]);
                }
            }


            // Wait for all to finish

            Pstream::waitRequests();

            // Collect neighbour fields

            for (label domain = 0; domain < Pstream::nProcs(); domain++)
            {
                const labelList& map = constructMap[domain];

                if (domain != Pstream::myProcNo() && map.size())
                {
                    if (recvFields[domain].size() != map.size())
                    {
                        FatalErrorIn
                        (
                            "template<class T>\n"
                            "void mapDistribute::distribute\n"
                            "(\n"
                            "    const Pstream::commsTypes commsType,\n"
                            "    const List<labelPair>& schedule,\n"
                            "    const label constructSize,\n"
                            "    const labelListList& subMap,\n"
                            "    const labelListList& constructMap,\n"
                            "    List<T>& field\n"
                            ")\n"
                        )   << "Expected from processor " << domain
                            << " " << map.size() << " but received "
                            << recvFields[domain].size() << " elements."
                            << abort(FatalError);
                    }

                    forAll(map, i)
                    {
                        cop(field[map[i]], recvFields[domain][i]);
                    }
                }
            }
        }
    }
    else
    {
        FatalErrorIn("mapDistribute::distribute(..)")
            << "Unknown communication schedule " << commsType
            << abort(FatalError);
    }
}


template<class T>
void Foam::mapDistribute::exchange
(
    const List<List<T> >& sendBuf,
    List<List<T> >& recvBuf
)
{
    if (!contiguous<T>())
    {
        FatalErrorIn("mapDistribute::exchange(..)")
            << "Not contiguous" << exit(FatalError);
    }

    if (Pstream::parRun())
    {

        // Determine sizes
        // ~~~~~~~~~~~~~~~

        labelListList allNTrans(Pstream::nProcs());
        allNTrans[Pstream::myProcNo()].setSize(Pstream::nProcs());

        forAll(allNTrans, procI)
        {
            allNTrans[Pstream::myProcNo()][procI] = sendBuf[procI].size();
        }
        combineReduce(allNTrans, listEq());


        // Set up receives
        // ~~~~~~~~~~~~~~~

        recvBuf.setSize(Pstream::nProcs());
        forAll(recvBuf, procI)
        {
            if (procI != Pstream::myProcNo())
            {
                recvBuf[procI].setSize(allNTrans[procI][Pstream::myProcNo()]);
                IPstream::read
                (
                    Pstream::nonBlocking,
                    procI,
                    reinterpret_cast<char*>(recvBuf[procI].begin()),
                    recvBuf[procI].byteSize()
                );
            }
        }

        // Set up sends
        // ~~~~~~~~~~~~

        forAll(sendBuf, procI)
        {
            if (procI != Pstream::myProcNo())
            {
                OPstream::write
                (
                    Pstream::nonBlocking,
                    procI,
                    reinterpret_cast<const char*>(sendBuf[procI].begin()),
                    sendBuf[procI].byteSize()
                );
            }
        }

        // Wait for completion
        Pstream::waitRequests();
    }

    // Do myself
    recvBuf[Pstream::myProcNo()] = sendBuf[Pstream::myProcNo()];
}


// ************************************************************************* //
