/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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

Application
    parallelTest

Description
    Test for various parallel routines.

\*---------------------------------------------------------------------------*/

#include "List.H"
#include "mapDistribute.H"
#include "argList.H"
#include "Time.H"
#include "IPstream.H"
#include "OPstream.H"
#include "vector.H"
#include "IOstreams.H"
#include "Random.H"
#include "Tuple2.H"

using namespace Foam;


void testMapDistribute()
{
    Random rndGen(43544*Pstream::myProcNo());

    // Generate random data.
    List<Tuple2<label, List<scalar>>> complexData(100);
    forAll(complexData, i)
    {
        complexData[i].first() = rndGen.position(0, Pstream::nProcs()-1);
        complexData[i].second().setSize(3);
        complexData[i].second()[0] = 1;
        complexData[i].second()[1] = 2;
        complexData[i].second()[2] = 3;
    }

    // Send all ones to processor indicated by .first()

    // Count how many to send
    labelList nSend(Pstream::nProcs(), 0);
    forAll(complexData, i)
    {
        label procI = complexData[i].first();
        nSend[procI]++;
    }

    // Collect items to be sent
    labelListList sendMap(Pstream::nProcs());
    forAll(sendMap, procI)
    {
        sendMap[procI].setSize(nSend[procI]);
    }
    nSend = 0;
    forAll(complexData, i)
    {
        label procI = complexData[i].first();
        sendMap[procI][nSend[procI]++] = i;
    }

    // Sync how many to send
    labelList nRecv;
    Pstream::exchangeSizes(sendMap, nRecv);

    // Collect items to be received
    labelListList recvMap(Pstream::nProcs());
    forAll(recvMap, procI)
    {
        recvMap[procI].setSize(nRecv[procI]);
    }

    label constructSize = 0;
    // Construct with my own elements first
    forAll(recvMap[Pstream::myProcNo()], i)
    {
        recvMap[Pstream::myProcNo()][i] = constructSize++;
    }
    // Construct from other processors
    forAll(recvMap, procI)
    {
        if (procI != Pstream::myProcNo())
        {
            forAll(recvMap[procI], i)
            {
                recvMap[procI][i] = constructSize++;
            }
        }
    }

    // Construct distribute map (destructively)
    mapDistribute map(constructSize, sendMap.xfer(), recvMap.xfer());

    // Distribute complexData
    map.distribute(complexData);

    Pout<< "complexData:" << complexData << endl;
}


// Print to Perr
template<class T>
Ostream& perrInfo(const T& data)
{
    Perr<< data;
    return Perr;
}


// Print to Perr
template<>
Ostream& perrInfo(const string& data)
{
    Perr<< data << " (size: " << data.size() << ")";
    return Perr;
}


template<class T>
void testTransfer(const T& input)
{
    T data = input;

    if (Pstream::master())
    {
        Perr<<"test transfer (" << (typeid(T).name()) << "): ";
        perrInfo(data) << nl << endl;
    }

    if (Pstream::myProcNo() != Pstream::masterNo())
    {
        {
            Perr<< "slave sending to master " << Pstream::masterNo() << endl;
            OPstream toMaster(Pstream::commsTypes::blocking, Pstream::masterNo());
            toMaster << data;
        }

        Perr<< "slave receiving from master " << Pstream::masterNo() << endl;
        IPstream fromMaster(Pstream::commsTypes::blocking, Pstream::masterNo());
        fromMaster >> data;
        perrInfo(data) << endl;
    }
    else
    {
        for
        (
            int slave = Pstream::firstSlave();
            slave <= Pstream::lastSlave();
            ++slave
        )
        {
            Perr<< "master receiving from slave " << slave << endl;
            IPstream fromSlave(Pstream::commsTypes::blocking, slave);
            fromSlave >> data;
            perrInfo(data) << endl;
        }

        for
        (
            int slave = Pstream::firstSlave();
            slave <= Pstream::lastSlave();
            ++slave
        )
        {
            Perr<< "master sending to slave " << slave << endl;
            OPstream toSlave(Pstream::commsTypes::blocking, slave);
            toSlave << data;
        }
    }
}


template<class T>
void testTokenized(const T& data)
{
    token tok;

    if (Pstream::master())
    {
        Perr<<"test tokenized \"" << data << "\"" << nl << endl;
    }

    if (Pstream::myProcNo() != Pstream::masterNo())
    {
        {
            Perr<< "slave sending to master " << Pstream::masterNo() << endl;
            OPstream toMaster
            (
                Pstream::commsTypes::blocking,
                Pstream::masterNo()
            );

            toMaster << data;
        }

        Perr<< "slave receiving from master " << Pstream::masterNo() << endl;
        IPstream fromMaster
        (
            Pstream::commsTypes::blocking,
            Pstream::masterNo()
        );

        fromMaster >> tok;
        Perr<< tok.info() << endl;
    }
    else
    {
        for
        (
            int slave = Pstream::firstSlave();
            slave <= Pstream::lastSlave();
            ++slave
        )
        {
            Perr<< "master receiving from slave " << slave << endl;
            IPstream fromSlave(Pstream::commsTypes::blocking, slave);
            fromSlave >> tok;
            Perr<< tok.info() << endl;
        }

        for
        (
            int slave = Pstream::firstSlave();
            slave <= Pstream::lastSlave();
            ++slave
        )
        {
            Perr<< "master sending to slave " << slave << endl;
            OPstream toSlave(Pstream::commsTypes::blocking, slave);
            toSlave << data;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"

    testMapDistribute();

    if (!Pstream::parRun())
    {
        Info<< "\nWarning: not parallel - skipping further tests\n" << endl;
        return 0;
    }

    Info<< "\nStarting transfers\n\n" << endl;

    testTransfer(vector(0, 1, 2));
    testTransfer(label(1234));
    testTransfer(scalar(3.14159));
    testTransfer(string("test   string"));
    testTransfer(string("  x "));

    {
        // Slightly roundabout way to construct with a nul in string
        string str1("embedded. nul character in string");
        str1[8] = '\0';

        Info<< "len: " << str1.size() << endl;
        testTransfer(str1);
    }
    testTransfer(word("3.141 59"));  // bad word, but transfer doesn't care

    testTokenized(label(1234));
    testTokenized(scalar(3.14159));
    testTokenized('a');
    testTokenized('$');  // will not tokenize well

    testTokenized(string("test   string1"));
    testTokenized("test   string1");
    testTokenized(word("3.141 59"));  // bad word, but transfer doesn't care

    testTokenized(string("  a "));
    testTokenized("  a ");

    testTokenized(string("  $ "));
    testTokenized("  $ ");  // reduces to 'char' and will not tokenize well

    testTokenized(string("  $$ "));
    testTokenized("  $$ "); // reduces to 'word' and is tagged as such


    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
