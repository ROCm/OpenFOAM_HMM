/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
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

Description

\*---------------------------------------------------------------------------*/

#include "IndirectList.H"
#include "IOstreams.H"
#include "Fstream.H"
#include "ListOps.H"
#include "labelIndList.H"
#include "argList.H"

using namespace Foam;

template<class ListType>
void printInfo(const ListType& lst)
{
    Info<< "full: " << flatOutput(lst.completeList()) << nl
        << "addr: " << flatOutput(lst.addressing()) << nl
        << "list: " << flatOutput(lst) << nl
        << endl;
}


template<class T, class ListType>
void testFind(const T& val, const ListType& lst)
{
    Info<< nl
        << "Search for "<< val << " in " << flatOutput(lst) << nl
        <<" found() = " << lst.found(val)
        <<" find() = " << lst.find(val)
        <<" rfind() = " << lst.rfind(val)
        <<" find(2) = " << lst.find(val, 2)
        <<" rfind(2) = " << lst.rfind(val, 2)
        <<" findIndex = " << findIndex(lst, val) << nl
        << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::addOption
    (
        "binary",
        "file",
        "write lists in binary to specified file"
    );

    argList args(argc, argv);

    List<label> completeList(20);

    forAll(completeList, i)
    {
        completeList[i] = 10*i;
    }

    Info<< "raw : " << flatOutput(completeList) << nl << endl;

    List<label> addresses{1, 0, 3, 7, 4, 8, 5, 1, 0, 3, 7, 4, 8, 5, };

    labelIndList idl1(completeList, addresses);

    printInfo(idl1);

    for (const label val : { 10, 30, 40, 50, 90, 80, 120 } )
    {
        testFind(val, idl1);
    }

    inplaceReverseList(addresses);

    idl1.resetAddressing(addresses.xfer());

    printInfo(idl1);

    // Test copying
    labelUIndList uidl1(idl1);
    labelIndList idl2(uidl1);
    labelIndList idl3(idl2);

    printInfo(uidl1);

    idl1.resetAddressing(List<label>());
//    idl2.resetAddressing(List<label>());

    Info<<"after resetAddressing:" << nl << endl;

    printInfo(uidl1);
    printInfo(idl1);
    printInfo(idl2);
    printInfo(idl3);

    fileName binaryOutput;
    if (args.optionReadIfPresent("binary", binaryOutput))
    {
        Info<<"Writing output to " << binaryOutput << endl;

        OFstream os(binaryOutput, IOstream::BINARY);

        os.writeEntry("idl1", idl1);
        os.writeEntry("idl2", idl2);
        os.writeEntry("idl3", idl3);
    }

    if (Pstream::parRun())
    {
        if (Pstream::master())
        {
            Pout<< "full: " << flatOutput(idl3.completeList()) << nl
                << "send: " << flatOutput(idl3) << endl;

            for
            (
                int slave = Pstream::firstSlave();
                slave <= Pstream::lastSlave();
                ++slave
            )
            {
                OPstream toSlave(Pstream::commsTypes::scheduled, slave);
                toSlave << idl3;
            }
        }
        else
        {
            // From master
            IPstream fromMaster
            (
                Pstream::commsTypes::scheduled,
                Pstream::masterNo()
            );

            List<label> recv(fromMaster);

            Pout<<"recv: " << flatOutput(recv) << endl;
        }

        // MPI barrier
        bool barrier = true;
        Pstream::scatter(barrier);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
