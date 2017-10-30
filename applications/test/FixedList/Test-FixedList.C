/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Test-FixedList

Description
    Simple tests and examples of use of FixedList

See also
    Foam::FixedList

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "FixedList.H"
#include "Fstream.H"
#include "List.H"
#include "IPstream.H"
#include "OPstream.H"
#include <numeric>

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList args(argc, argv);

    if (false)
    {
        FixedList<string, 1> ident;

        auto iter = ident.begin();

        Info << iter->size() << endl;

        auto riter = ident.rbegin();
        Info << riter->size() << endl;

        auto iter2 = ident.rbegin();

        iter2 = iter;
    }

    {
        FixedList<label, 15> ident;
        std::iota(ident.begin(), ident.end(), 0);

        // auto iter = ident.begin();
        //
        // iter += 5;
        // Info << *iter << "< " << endl;
        // iter -= 2;
        // Info << *iter << "< " << endl;

        // Don't yet bother with making reverse iterators random access
        // auto riter = ident.crbegin();

        // riter += 5;
        // Info << *riter << "< " << endl;
        // riter += 2;
        // Info << *riter << "< " << endl;

        Info<<"Ident:";
        forAllConstIters(ident, iter)
        {
            Info<<" " << *iter;
        }
        Info<< nl;

        Info<<"reverse:";
        forAllReverseIters(ident, iter)
        {
            Info<<" " << *iter;
        }
        Info<< nl;

        Info<<"const reverse:";
        forAllConstReverseIters(ident, iter)
        {
            Info<<" " << *iter;
        }
        Info<< nl;
    }

    {
        FixedList<label, 4> list1{1, 2, 3, 4};

        Info<< "list1:" << list1
            << " hash:" << FixedList<label, 4>::Hash<>()(list1) << endl;

        label a[4] = {0, 1, 2, 3};
        FixedList<label, 4> list2(a);

        Info<< "list2:" << list2
            << " hash:" << FixedList<label, 4>::Hash<>()(list2) << endl;

        // Using FixedList for content too
        {
            List<FixedList<label, 4>> twolists{list1, list2};
            Info<<"List of FixedList: " << flatOutput(twolists) << endl;
            sort(twolists);
            // outer-sort only
            Info<<"sorted FixedList : " << flatOutput(twolists) << endl;
        }

        Info<< "====" << nl
            << "Test swap" << nl;

        Info<< "list1: " << list1 << nl
            << "list2: " << list2 << endl;
        list1.swap(list2);
        Info<< "The swap() method" << endl;
        Info<< "list1: " << list1 << nl
            << "list2: " << list2 << endl;

        Swap(list1, list2);
        Info<< "The Swap() function" << endl;
        Info<< "list1: " << list1 << nl
            << "list2: " << list2 << endl;

        Info<< "====" << nl;
    }

    List<label> list3{0, 1, 2, 3};
    FixedList<label, 4> list4(list3.begin(), list3.end());
    Info<< "list3: " << list3 << nl
        << "list4: " << list4 << endl;

    list4 = {1, 2, 3, 5};
    Info<< "list4: " << list4 << nl;

    FixedList<label, 5> list5{0, 1, 2, 3, 4};
    Info<< "list5: " << list5 << endl;

    List<FixedList<label, 2>> list6{{0, 1}, {2, 3}};
    Info<< "list6: " << list6 << endl;

    if (Pstream::parRun())
    {
        if (Pstream::myProcNo() != Pstream::masterNo())
        {
            Serr<< "slave sending to master "
                << Pstream::masterNo() << endl;

            OPstream toMaster
            (
                Pstream::commsTypes::blocking,
                Pstream::masterNo()
            );

            FixedList<label, 2> list3;
            list3[0] = 0;
            list3[1] = 1;
            toMaster << list3;
        }
        else
        {
            for
            (
                int slave = Pstream::firstSlave();
                slave <= Pstream::lastSlave();
                slave++
            )
            {
                Serr << "master receiving from slave " << slave << endl;
                IPstream fromSlave(Pstream::commsTypes::blocking, slave);
                FixedList<label, 2> list3(fromSlave);

                Serr<< list3 << endl;
            }
        }
    }

    return 0;
}


// ************************************************************************* //
