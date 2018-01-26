/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

Description

\*---------------------------------------------------------------------------*/

#include "OSspecific.H"

#include "IOstreams.H"
#include "labelList.H"
#include "DynamicList.H"
#include "face.H"
#include "pointField.H"
#include "DynamicField.H"

using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    labelList lstA(10);
    labelList lstC
    {
        1, 2, 3, 4
    };

    forAll(lstA, i)
    {
        lstA[i] = 5 - i;
    }

    Info<< "lstA: " << lstA << nl
        << "lstC: " << lstC << nl;

    Xfer<labelList> xA = xferMove(lstA);
    Xfer<labelList> xB;

    labelList lstB( xA );

    Info<< "xA: " << xA() << nl
        << "xB: " << xB() << nl
        << "lstA: " << lstA << nl
        << "lstB: " << lstB << nl
        << "lstC: " << lstC << nl;

    // Now illegal: xA = lstB;

    xA->transfer(lstB);

    Info<< "xA: " << xA() << nl
        << "xB: " << xB() << nl
        << "lstA: " << lstA << nl
        << "lstB: " << lstB << nl
        << "lstC: " << lstC << nl;

    xB = xA;

    // Construct with forwarding. For this example, truly ugly.
    Xfer<labelList> xFwdA =
        Xfer<labelList>::New
        (
            std::initializer_list<label>
            {
               1, 2, 10, 20, 15, 24, 200
            }
        );

    Xfer<labelList> xFwdB = Xfer<labelList>::New(label(8), 123);
    Xfer<labelList> xFwdC = Xfer<labelList>::New();

    Info<< nl
        << "Constructed with forwarding: " << nl
        << *xFwdA << nl
        << *xFwdB << nl
        << *xFwdC << nl
        << nl;


    labelList lstD(xferCopy(lstC));
    labelList lstE(xferMove(lstC));

    // this must be empty
    labelList lstF = xferCopy(lstC);

    Info<< "xA: " << xA() << nl
        << "xB: " << xB() << nl
        << "lstA: " << lstA << nl
        << "lstB: " << lstB << nl
        << "lstC: " << lstC << nl
        << "lstD: " << lstD << nl
        << "lstE: " << lstE << nl
        << "lstF: " << lstF << nl;

    Info<< "xB[" << xB->size() << "]\n";

    // clear the underlying List
    xB->clear();

    Info<< "xB[" << xB->size() << "]\n";

    DynamicList<label> dl;
    for (label i = 0; i < 5; ++i)
    {
        dl.append(i);
    }

    face f1(dl);
    face f2(xferCopy<labelList>(dl));

    Info<< "dl[" << dl.size() << "/" << dl.capacity() << "] " << dl << nl;
    Info<< "f1: " << f1 << nl;
    Info<< "f2: " << f2 << nl;

    // add some more labels
    for (label i = 5; i < 8; ++i)
    {
        dl.append(i);
    }

    // note: xfer() method returns a plain labelList
    face f3(dl.xfer());
    Info<< "dl[" << dl.size() << "/" << dl.capacity() << "] " << dl << nl;
    Info<< "f3: " << f3 << nl;

    Info<<"\nflip faces:" << nl;
    f1.flip();
    f3.flip();
    Info<< "f1: " << f1 << nl;
    Info<< "f3: " << f3 << nl;


    {
        Info<<"\nTest xfer with fields:" << nl;
        List<point> list1
        {
            { 0, 1, 2 },
            { 3, 4, 5 },
            { 6, 7, 8 },
            { 9, 10, 11 },
        };

        // Field from Xfer<List>
        pointField field1(list1.xfer());
        Info<<nl
            << "xfer construct from List" << nl
            <<"input (list) = " << list1 << nl
            <<"output (field) = " << field1 << nl;


        // Field from Xfer<List> ... again
        pointField field2(field1.xfer());
        Info<<nl
            <<"xfer construct from Field (as List): " << nl
            <<"input (field) = " << field1 << nl
            <<"output (field) = " << field2 << nl;


        // Field from Xfer<Field>
        pointField field3(xferMove(field2));
        Info<<nl
            <<"xfer construct from Field (as Field): " << nl
            <<"input (field) = " << field2 << nl
            <<"output (field) = " << field3 << nl;


        // Field from Xfer<Field> .. again
        pointField field4(xferCopy(field3));
        Info<<nl
            <<"xfer copy construct from Field (as Field): " << nl
            <<"input (field) = " << field3 << nl
            <<"output (field) = " << field4 << nl;


        DynamicField<point> dyfield1(xferCopy(field4));
        Info<<nl
            <<"xfer copy construct from Field (as Field): " << nl
            <<"input (field) = " << field4 << nl
            <<"output (dyn-field) = " << dyfield1 << nl;
    }

    return 0;
}


// ************************************************************************* //
