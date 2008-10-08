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
    Testing TIP4P potential and dynamics

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
//#include "md.H"
#include "diagTensor.H"

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

    qSites = qSites*1.602176487e-19;

    scalar rH = 0.9572e-10;
    scalar rM = 0.15e-10;
    scalar thetaH = mathematicalConstant::pi*104.52/(2.0*180.0);

    List<vector> pSites(nSites);

    pSites[0] = vector
    (
        rH*Foam::sin(thetaH),
        rH*Foam::cos(thetaH)*(1.0 - 1.0/(1.0 + 0.5*mSites[2]/mSites[0])),
        0
    );
    pSites[1] = vector(-pSites[0].x(), pSites[0].y(), 0);
    pSites[2] = vector(0, -Foam::cos(thetaH)*rH/(1.0 +  0.5*mSites[2]/mSites[0]), 0);
    pSites[3] = vector(0, rM + pSites[2].y(), 0);

    diagTensor I
    (
        mSites[2] * pSites[2].y() * pSites[2].y()
          + 2.0 * mSites[0] * pSites[0].y() * pSites[0].y(),
        2.0 * mSites[0] + pSites[0].x() * pSites[0].x(),
        mSites[2] * pSites[2].y() * pSites[2].y()
          + 2.0 * mSites[0] * (pSites[0].y() * pSites[0].y() + pSites[0].x() * pSites[0].x())
    );

    Info<< m
        << nl << qSites
        << nl << pSites
        << nl << I << endl;

    vector p1(0, 0, 0);

    vector v1(100, 0, 0);

    tensor R1(tensor(1, 0, 0, 0, 1, 0, 0, 0, 1));

    vector omega1(0, -2e12, 0);

    List<vector> pSites1 = R1 & pSites;

    List<vector> fSites1(nSites, vector::zero);

    vector a1(vector::zero);

    vector alpha1(vector::zero);

    Info<< nl << p1
        << nl << v1
        << nl << R1
        << nl << omega1
        << nl << fSites1
        << nl << pSites1
        << endl;

    // scalar freq = 1e12;
    // scalar ampl = 3e7;

    vector eForce(vector::zero);

    Info<< "\nStarting time loop\n" << endl;

    Info<< pSites.size() << nl << "Water test"
        << nl << "H1"
        << " " << pSites1[0].x() << " " << pSites1[0].y() << " " << pSites1[0].z()
        << nl << "H2"
        << " " << pSites1[1].x() << " " << pSites1[1].y() << " " << pSites1[1].z()
        << nl << "0"
        << " " << pSites1[2].x() << " " << pSites1[2].y() << " " << pSites1[2].z()
        << nl << "M"
        << " " << pSites1[3].x() << " " << pSites1[3].y() << " " << pSites1[3].z()
        << endl;

    while (runTime.run())
    {
        runTime++;

        a1 = vector(-5e13, 0, 0);
        alpha1 = vector(0, 1e24, 0);

        // Leapfrog part 1
        v1 += 0.5*runTime.deltaT().value()*a1;

        p1 += runTime.deltaT().value()*v1;

        omega1 += 0.5*runTime.deltaT().value()*alpha1;

        scalar phi1 = 0.5*runTime.deltaT().value()*omega1.x();

        tensor U1
        (
            tensor
            (
                1, 0, 0,
                0, Foam::cos(phi1), Foam::sin(phi1),
                0, -Foam::sin(phi1), Foam::cos(phi1)
            )
        );

        scalar phi2 = 0.5*runTime.deltaT().value()*omega1.y();

        tensor U2
        (
            tensor
            (
                Foam::cos(phi2), 0, -Foam::sin(phi2),
                0, 1, 0,
                Foam::sin(phi2), 0, Foam::cos(phi2)
            )
        );

        scalar phi3 = runTime.deltaT().value()*omega1.z();

        tensor U3
        (
            tensor
            (
                Foam::cos(phi3), Foam::sin(phi3), 0,
                -Foam::sin(phi3), Foam::cos(phi3), 0,
                0, 0, 1
            )
        );

        tensor UT = U1.T() & U2.T() & U3.T() & U2.T() & U1.T();

        R1 = R1 & UT;

        pSites1 = p1 + (R1 & pSites);

        // Force calculation

        a1 = vector();

        //        eForce.z() = ampl*Foam::sin(2*mathematicalConstant::pi*freq*runTime.timeOutputValue());

        //        fSites1 = qSites * eForce;

        vector tau1(vector::zero);

        forAll(fSites1, fS)
        {
            const vector& f = fSites1[fS];

            a1 += f/m;

            tau1 += (pSites1[fS] - p1) ^ f;
        }

        alpha1 = R1 & inv(I) & R1.T() & tau1;

        a1 = vector(-5e13, 0, 0);
        alpha1 = vector(0, 1e24, 0);

        // Leapfrog part 2

        v1 += 0.5*runTime.deltaT().value()*a1;

        omega1 += 0.5*runTime.deltaT().value()*alpha1;

        Info<< pSites.size() << nl << "Water test"
            << nl << "H1"
            << " " << pSites1[0].x() << " " << pSites1[0].y() << " " << pSites1[0].z()
            << nl << "H2"
            << " " << pSites1[1].x() << " " << pSites1[1].y() << " " << pSites1[1].z()
            << nl << "0"
            << " " << pSites1[2].x() << " " << pSites1[2].y() << " " << pSites1[2].z()
            << nl << "M"
            << " " << pSites1[3].x() << " " << pSites1[3].y() << " " << pSites1[3].z()
            << endl;

        // Info<< "Time = " << runTime.timeName() << endl;
    }

    Info<< v1 << nl << omega1 << endl;

    Info << "End\n" << endl;

    return(0);
}

