/*---------------------------------------------------------------------------*\
 =========                   |
 \\      /   F ield          | OpenFOAM: The Open Source CFD Toolbox
  \\    /    O peration      |
   \\  /     A nd            | Copyright (C) 2008-2009 OpenCFD Ltd.
    \\/      M anipulation   |
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

Application
    DistributionTest

Description
    Test the Distribution class

\*---------------------------------------------------------------------------*/

#include "vector.H"
#include "labelVector.H"
#include "tensor.H"
#include "Distribution.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

int main(int argc, char *argv[])
{
    Distribution<scalar> dS(scalar(3.0));
    Distribution<vector> dV(vector(0.1, 0.24, 0.18));
    Distribution<labelVector> dLV(labelVector(2,3,4));

    dS.add(1.0);
    dS.add(2.0);
    dS.add(12.0);
    dS.add(1.3);

    Info<< dS << nl << dS.raw() << endl;

    vector vA(1.2, 1.3, 1.1);
    vector vB(1.3, 1.5, 1.6);
    vector vC(0.5, 5.3, 1.1);

    dV.add(vA);
    dV.add(vB);
    dV.add(vC);

    Info<< dV << nl << dV.raw() << endl;

    labelVector lVA(6, 8, 11);
    labelVector lVB(-12, -3, 6);
    labelVector lVC(-4, -2, 5);

    dLV.add(lVA);
    dLV.add(lVB);
    dLV.add(lVC);

    Info<< dLV << nl << dLV.raw() << endl;

    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
