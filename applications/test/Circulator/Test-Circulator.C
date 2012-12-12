/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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
    Test-circulator

Description

\*---------------------------------------------------------------------------*/

#include "List.H"
#include "ListOps.H"
#include "face.H"
#include "circulator.H"
#include "const_circulator.H"


using namespace Foam;


// return
//   0: no match
//  +1: identical
//  -1: same face, but different orientation
label compare(const face& a, const face& b)
{
    // Basic rule: we assume that the sequence of labels in each list
    // will be circular in the same order (but not necessarily in the
    // same direction or from the same starting point).

    // Trivial reject: faces are different size
    label sizeA = a.size();
    label sizeB = b.size();

    if (sizeA != sizeB || sizeA == 0)
    {
        return 0;
    }

    const_circulator<face> aCirc(a);
    const_circulator<face> bCirc(b);

    // Rotate face b until its element matches the starting element of face a.
    do
    {
        if (aCirc() == bCirc())
        {
            // Set bCirc fulcrum to its iterator and increment the iterators
            bCirc.setFulcrumToIterator();
            ++aCirc;
            ++bCirc;

            break;
        }
    } while (bCirc.circulate(CirculatorBase::CLOCKWISE));

    // Look forwards around the faces for a match
    do
    {
        if (aCirc() != bCirc())
        {
            break;
        }
    }
    while
    (
        aCirc.circulate(CirculatorBase::CLOCKWISE),
        bCirc.circulate(CirculatorBase::CLOCKWISE)
    );

    // If the circulator has stopped then faces a and b matched.
    if (!aCirc.circulate())
    {
        return 1;
    }
    else
    {
        // Reset the circulators back to their fulcrum
        aCirc.setIteratorToFulcrum();
        bCirc.setIteratorToFulcrum();
        ++aCirc;
        --bCirc;
    }

    // Look backwards around the faces for a match
    do
    {
        if (aCirc() != bCirc())
        {
            break;
        }
    }
    while
    (
        aCirc.circulate(CirculatorBase::CLOCKWISE),
        bCirc.circulate(CirculatorBase::ANTICLOCKWISE)
    );

    // If the circulator has stopped then faces a and b matched.
    if (!aCirc.circulate())
    {
        return -1;
    }

    return 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    Info<< "Test the implementation of a circular iterator" << nl << endl;

    Info<< "Test const circulator. First go forwards, then backwards."
        << nl << endl;

    face f(identity(4));

    const_circulator<face> cStart(f);

    if (cStart.size()) do
    {
        Info<< "Iterate forwards over face (prev/curr/next) : "
            << cStart.prev() << " / " << cStart() << " / " << cStart.next()
            << endl;

    } while (cStart.circulate(CirculatorBase::CLOCKWISE));

    if (cStart.size()) do
    {
        Info<< "Iterate backwards over face : " << cStart() << endl;

    } while (cStart.circulate(CirculatorBase::ANTICLOCKWISE));


    Info<< nl << nl << "Test non-const circulator" << nl << endl;

    circulator<face> cStart2(f);

    Info<< "Face before : " << f << endl;

    if (cStart2.size()) do
    {
        Info<< "Iterate forwards over face (prev/curr/next) : "
            << cStart2.prev() << " / " << cStart2() << " / " << cStart2.next()
            << endl;

    } while (cStart2.circulate(CirculatorBase::CLOCKWISE));

    if (cStart2.size()) do
    {
        Info<< "Iterate forwards over face, adding 1 to each element : "
            << cStart2();

        cStart2() += 1;

        Info<< " -> " << cStart2() << endl;
    } while (cStart2.circulate(CirculatorBase::CLOCKWISE));

    Info<< "Face after : " << f << endl;


    Info<< nl << nl << "Compare two faces: " << endl;
    face a(identity(5));
    Info<< "Compare " << a << " and " << a << " Match = " << compare(a, a)
        << endl;

    face b(reverseList(a));
    Info<< "Compare " << a << " and " << b << " Match = " << compare(a, b)
        << endl;

    face c(a);
    c[4] = 3;
    Info<< "Compare " << a << " and " << c << " Match = " << compare(a, c)
        << endl;

    face d(rotateList(a, 2));
    Info<< "Compare " << a << " and " << d << " Match = " << compare(a, d)
        << endl;

    face g(labelList(5, 1));
    face h(g);
    Info<< "Compare " << g << " and " << h << " Match = " << compare(g, h)
        << endl;

    g[0] = 2;
    h[3] = 2;
    Info<< "Compare " << g << " and " << h << " Match = " << compare(g, h)
        << endl;

    g[4] = 3;
    h[4] = 3;
    Info<< "Compare " << g << " and " << h << " Match = " << compare(g, h)
        << endl;

    face face1(identity(1));
    Info<< "Compare " << face1 << " and " << face1
        << " Match = " << compare(face1, face1) << endl;

    Info<< nl << nl << "Zero face" << nl << endl;

    face fZero;
    circulator<face> cZero(fZero);

    if (cZero.size()) do
    {
        Info<< "Iterate forwards over face : " << cZero() << endl;

    } while (cZero.circulate(CirculatorBase::CLOCKWISE));

    fZero = face(identity(5));

    // circulator was invalidated so reset
    cZero = circulator<face>(fZero);

    do
    {
        Info<< "Iterate forwards over face : " << cZero() << endl;

    } while (cZero.circulate(CirculatorBase::CLOCKWISE));


    Info<< nl << nl << "Simultaneously go forwards/backwards over face " << f
        << nl << endl;

    const_circulator<face> circForward(f);
    const_circulator<face> circBackward(f);

    if (circForward.size() && circBackward.size()) do
    {
        Info<< "Iterate over face forwards : " << circForward()
            << ", backwards : " << circBackward() << endl;
    }
    while
    (
        circForward.circulate(CirculatorBase::CLOCKWISE),
        circBackward.circulate(CirculatorBase::ANTICLOCKWISE)
    );

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
