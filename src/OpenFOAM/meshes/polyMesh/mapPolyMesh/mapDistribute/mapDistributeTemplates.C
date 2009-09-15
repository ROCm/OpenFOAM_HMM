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
            // 1. convert to contiguous buffer
            // 2. send buffer
            // 3. receive buffer
            // 4. read from buffer into List<T>

            List<List<char> > sendFields(Pstream::nProcs());
            labelListList allNTrans(Pstream::nProcs());
            labelList& nsTransPs = allNTrans[Pstream::myProcNo()];
            nsTransPs.setSize(Pstream::nProcs(), 0);

            // Stream data into sendField buffers
            for (label domain = 0; domain < Pstream::nProcs(); domain++)
            {
                const labelList& map = subMap[domain];

                if (domain != Pstream::myProcNo() && map.size())
                {
                    // Put data into send buffer
                    OPstream toDomain(Pstream::nonBlocking, domain);
                    toDomain << UIndirectList<T>(field, map);

                    // Store the size
                    nsTransPs[domain] = toDomain.bufPosition();

                    // Transfer buffer out
                    sendFields[domain].transfer(toDomain.buf());
                    toDomain.bufPosition() = 0;

                }
            }

            // Send sizes across
            combineReduce(allNTrans, listEq());

            // Start sending buffers
            for (label domain = 0; domain < Pstream::nProcs(); domain++)
            {
                const labelList& map = subMap[domain];

                if (domain != Pstream::myProcNo() && map.size())
                {
                    OPstream::write
                    (
                        Pstream::nonBlocking,
                        domain,
                        reinterpret_cast<const char*>
                        (
                            sendFields[domain].begin()
                        ),
                        nsTransPs[domain]
                    );
                }
            }

            // Set up receives from neighbours

            PtrList<IPstream> fromSlave(Pstream::nProcs());

            for (label domain = 0; domain < Pstream::nProcs(); domain++)
            {
                const labelList& map = constructMap[domain];

                if (domain != Pstream::myProcNo() && map.size())
                {
                    // Start receiving
                    fromSlave.set
                    (
                        domain,
                        new IPstream
                        (
                            Pstream::nonBlocking,
                            domain,
                            allNTrans[domain][Pstream::myProcNo()]
                        )
                    );
                }
            }


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


            // Wait till all finished
            Pstream::waitRequests();

            // Consume
            for (label domain = 0; domain < Pstream::nProcs(); domain++)
            {
                const labelList& map = constructMap[domain];

                if (domain != Pstream::myProcNo() && map.size())
                {
                    List<T> recvField(fromSlave[domain]);

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

                    // Delete receive buffer
                    fromSlave.set(domain, NULL);
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
            // 1. convert to contiguous buffer
            // 2. send buffer
            // 3. receive buffer
            // 4. read from buffer into List<T>

            List<List<char> > sendFields(Pstream::nProcs());
            labelListList allNTrans(Pstream::nProcs());
            labelList& nsTransPs = allNTrans[Pstream::myProcNo()];
            nsTransPs.setSize(Pstream::nProcs());

            // Stream data into sendField buffers
            for (label domain = 0; domain < Pstream::nProcs(); domain++)
            {
                const labelList& map = subMap[domain];

                if (domain != Pstream::myProcNo() && map.size())
                {
                    // Put data into send buffer
                    OPstream toDomain(Pstream::nonBlocking, domain);
                    toDomain << UIndirectList<T>(field, map);

                    // Store the size
                    nsTransPs[domain] = toDomain.bufPosition();

                    // Transfer buffer out
                    sendFields[domain].transfer(toDomain.buf());
                    toDomain.bufPosition() = 0;
                }
            }

            // Send sizes across
            combineReduce(allNTrans, listEq());

            // Start sending buffers
            for (label domain = 0; domain < Pstream::nProcs(); domain++)
            {
                const labelList& map = subMap[domain];

                if (domain != Pstream::myProcNo() && map.size())
                {
                    OPstream::write
                    (
                        Pstream::nonBlocking,
                        domain,
                        reinterpret_cast<const char*>
                        (
                            sendFields[domain].begin()
                        ),
                        nsTransPs[domain]
                    );
                }
            }

            // Set up receives from neighbours

            PtrList<IPstream> fromSlave(Pstream::nProcs());

            for (label domain = 0; domain < Pstream::nProcs(); domain++)
            {
                const labelList& map = constructMap[domain];

                if (domain != Pstream::myProcNo() && map.size())
                {
                    // Start receiving
                    fromSlave.set
                    (
                        domain,
                        new IPstream
                        (
                            Pstream::nonBlocking,
                            domain,
                            allNTrans[domain][Pstream::myProcNo()]
                        )
                    );
                }
            }


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
            Pstream::waitRequests();

            // Consume
            for (label domain = 0; domain < Pstream::nProcs(); domain++)
            {
                const labelList& map = constructMap[domain];

                if (domain != Pstream::myProcNo() && map.size())
                {
                    List<T> recvField(fromSlave[domain]);

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

                    // Delete receive buffer
                    fromSlave.set(domain, NULL);
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
                    IPstream::read
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


// ************************************************************************* //
