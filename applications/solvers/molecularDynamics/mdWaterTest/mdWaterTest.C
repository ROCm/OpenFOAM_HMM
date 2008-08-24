/*---------------------------------------------------------------------------*\
 =========                   |
 \\      /   F ield          | OpenFOAM: The Open Source CFD Toolbox
  \\    /    O peration      |
   \\  /     A nd            | Copyright (C) 1991-2008 OpenCFD Ltd.
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
    mdWaterTest
Description
    Testing TIP4P potential and dynamics test

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "md.H"

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    label nSites = 4;

    scalarList mSites(nSites);

    mSites[0] = 1.67353255e-27;
    mSites[1] = mSites[0];
    mSites[2] = 2.6560176e-26;
    mSites[3] = 0.0;

    scalar m = sum(mSites);

    scalarList qSites(nSites);

    qSites[0] = 0.52;
    qSites[1] = qSites[0];
    qSites[2] = 0.0;
    qSites[3] = -2.0*qSites[0];

    qSites *= 1.602176487e-19;

    scalar rH = 0.9572e-10;
    scalar rM = 0.15e-10;
    scalar thetaH = 104.52/2.0;

    vectorList pSites(nSites);

    pSites[0] = rH*vector
    (
        sin(thetaH),
        cos(thetaH)*(1.0 - 1.0/(1.0 + 0.5*mSites[2]/mSites[0))),
        0
    );
    pSites[1] = vector(-pSites[0].x(), pSites[0].x(), 0);
    pSites[2] = vector(0, -cos(thetaH)*rH/(1.0 +  0.5*mSites[2]/mSites[0]), 0);
    ps[3] = vector(0, rM + pSites[2].y(), 0);

    diagTensor I
    (
        mSites[2] * pSites[2].y() * pSites[2].y()
          + 2.0 * mSites[0] * pSites[0].y() * pSites[0].y(),
        2.0 * mSites[0] + pSites[0].x() * pSites[0].x(),
        mSites[2] * pSites[2].y() * pSites[2].y()
          + 2.0 * mSites[0] * (pSites[0].y() * pSites[0].y() + pSites[0].x() * pSites[0].x())
    );

    vector p1(0, 0, 0);

    vector v1(150.0, -23.0, 92);

    tensor R1(tensor::one);

    vector omega1(1.2e13, -2.6e11, 5.9e12);

    vectorList fSites1(nSites, vector::zero);

    Info << "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        runTime++;

        Info << "Time = " << runTime.timeName() << endl;
    }

    Info << "End\n" << endl;

    return(0);
}

