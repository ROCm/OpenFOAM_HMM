/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2007-2009 OpenCFD Ltd.
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

Class
    momentOfInertia

Description
    Reimplementation of volInt.c by Brian Mirtich.
        *  mirtich@cs.berkeley.edu                             *
        *  http://www.cs.berkeley.edu/~mirtich                 *

-------------------------------------------------------------------------------
*/

#include "momentOfInertia.H"
//#include "pyramidPointFaceRef.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

//Foam::tensor Foam::momentOfInertia
//(
//    const pointField& points,
//    const faceList& faces,
//    const cell& cFaces,
//    const point& cc
//)
//{
//    tensor t(tensor::zero);
//
//    forAll(cFaces, i)
//    {
//        const face& f = faces[cFaces[i]];
//
//        scalar pyrVol = pyramidPointFaceRef(f, cc).mag(points);
//
//        vector pyrCentre = pyramidPointFaceRef(f, cc).centre(points);
//
//        vector d = pyrCentre - cc;
//
//        t.xx() += pyrVol*(sqr(d.y()) + sqr(d.z()));
//        t.yy() += pyrVol*(sqr(d.x()) + sqr(d.z()));
//        t.zz() += pyrVol*(sqr(d.x()) + sqr(d.y()));
//
//        t.xy() -= pyrVol*d.x()*d.y();
//        t.xz() -= pyrVol*d.x()*d.z();
//        t.yz() -= pyrVol*d.y()*d.z();
//    }
//
//    // Symmetric
//    t.yx() = t.xy();
//    t.zx() = t.xz();
//    t.zy() = t.yz();
//
//    return t;
//}


#define sqr(x) ((x)*(x))
#define pow3(x) ((x)*(x)*(x))

// compute various integrations over projection of face
void Foam::compProjectionIntegrals
(
    const pointField& points,
    const face& f,
    const direction A,
    const direction B,

    scalar& P1,
    scalar& Pa,
    scalar& Pb,
    scalar& Paa,
    scalar& Pab,
    scalar& Pbb,
    scalar& Paaa,
    scalar& Paab,
    scalar& Pabb,
    scalar& Pbbb
)
{
    P1 = Pa = Pb = Paa = Pab = Pbb = Paaa = Paab = Pabb = Pbbb = 0.0;

    forAll(f, i)
    {
        scalar a0 = points[f[i]][A];
        scalar b0 = points[f[i]][B];
        scalar a1 = points[f[(i+1) % f.size()]][A];
        scalar b1 = points[f[(i+1) % f.size()]][B];
        scalar da = a1 - a0;
        scalar db = b1 - b0;

        scalar a0_2 = a0 * a0;
        scalar a0_3 = a0_2 * a0;
        scalar a0_4 = a0_3 * a0;

        scalar b0_2 = b0 * b0;
        scalar b0_3 = b0_2 * b0;
        scalar b0_4 = b0_3 * b0;

        scalar a1_2 = a1 * a1;
        scalar a1_3 = a1_2 * a1;

        scalar b1_2 = b1 * b1;
        scalar b1_3 = b1_2 * b1;

        scalar C1 = a1 + a0;

        scalar Ca = a1*C1 + a0_2;
        scalar Caa = a1*Ca + a0_3;
        scalar Caaa = a1*Caa + a0_4;

        scalar Cb = b1*(b1 + b0) + b0_2;
        scalar Cbb = b1*Cb + b0_3;
        scalar Cbbb = b1*Cbb + b0_4;

        scalar Cab = 3*a1_2 + 2*a1*a0 + a0_2;
        scalar Kab = a1_2 + 2*a1*a0 + 3*a0_2;

        scalar Caab = a0*Cab + 4*a1_3;
        scalar Kaab = a1*Kab + 4*a0_3;

        scalar Cabb = 4*b1_3 + 3*b1_2*b0 + 2*b1*b0_2 + b0_3;
        scalar Kabb = b1_3 + 2*b1_2*b0 + 3*b1*b0_2 + 4*b0_3;

        P1 += db*C1;
        Pa += db*Ca;
        Paa += db*Caa;
        Paaa += db*Caaa;
        Pb += da*Cb;
        Pbb += da*Cbb;
        Pbbb += da*Cbbb;
        Pab += db*(b1*Cab + b0*Kab);
        Paab += db*(b1*Caab + b0*Kaab);
        Pabb += da*(a1*Cabb + a0*Kabb);
  }

  P1 /= 2.0;
  Pa /= 6.0;
  Paa /= 12.0;
  Paaa /= 20.0;
  Pb /= -6.0;
  Pbb /= -12.0;
  Pbbb /= -20.0;
  Pab /= 24.0;
  Paab /= 60.0;
  Pabb /= -60.0;
}


void Foam::compFaceIntegrals
(
    const pointField& points,
    const face& f,
    const vector& n,
    const scalar w,
    const direction A,
    const direction B,
    const direction C,

    scalar& Fa,
    scalar& Fb,
    scalar& Fc,
    scalar& Faa,
    scalar& Fbb,
    scalar& Fcc,
    scalar& Faaa,
    scalar& Fbbb,
    scalar& Fccc,
    scalar& Faab,
    scalar& Fbbc,
    scalar& Fcca
)
{
    scalar P1, Pa, Pb, Paa, Pab, Pbb, Paaa, Paab, Pabb, Pbbb;

    compProjectionIntegrals
    (
        points,
        f,
        A,
        B,

        P1,
        Pa,
        Pb,
        Paa,
        Pab,
        Pbb,
        Paaa,
        Paab,
        Pabb,
        Pbbb
    );

    scalar k1 = 1 / n[C];
    scalar k2 = k1 * k1;
    scalar k3 = k2 * k1;
    scalar k4 = k3 * k1;

    Fa = k1 * Pa;
    Fb = k1 * Pb;
    Fc = -k2 * (n[A]*Pa + n[B]*Pb + w*P1);

    Faa = k1 * Paa;
    Fbb = k1 * Pbb;
    Fcc = k3 * (sqr(n[A])*Paa + 2*n[A]*n[B]*Pab + sqr(n[B])*Pbb
         + w*(2*(n[A]*Pa + n[B]*Pb) + w*P1));

    Faaa = k1 * Paaa;
    Fbbb = k1 * Pbbb;
    Fccc = -k4 * (pow3(n[A])*Paaa + 3*sqr(n[A])*n[B]*Paab
           + 3*n[A]*sqr(n[B])*Pabb + pow3(n[B])*Pbbb
           + 3*w*(sqr(n[A])*Paa + 2*n[A]*n[B]*Pab + sqr(n[B])*Pbb)
           + w*w*(3*(n[A]*Pa + n[B]*Pb) + w*P1));

    Faab = k1 * Paab;
    Fbbc = -k2 * (n[A]*Pabb + n[B]*Pbbb + w*Pbb);
    Fcca = k3 * (sqr(n[A])*Paaa + 2*n[A]*n[B]*Paab + sqr(n[B])*Pabb
         + w*(2*(n[A]*Paa + n[B]*Pab) + w*Pa));
}


void Foam::compVolumeIntegrals
(
    const pointField& points,
    const faceList& faces,
    const cell& cFaces,
    const vectorField& fNorm,
    const scalarField& fW,

    scalar& T0,
    vector& T1,
    vector& T2,
    vector& TP
)
{
    T0 = 0;
    T1 = vector::zero;
    T2 = vector::zero;
    TP = vector::zero;

    forAll(cFaces, i)
    {
        const vector& n = fNorm[i];

        scalar nx = mag(n[0]);
        scalar ny = mag(n[1]);
        scalar nz = mag(n[2]);

        direction A, B, C;

        if (nx > ny && nx > nz)
        {
            C = 0;
        }
        else
        {
            C = (ny > nz) ? 1 : 2;
        }

        A = (C + 1) % 3;
        B = (A + 1) % 3;

        scalar Fa, Fb, Fc, Faa, Fbb, Fcc, Faaa, Fbbb, Fccc, Faab, Fbbc, Fcca;
        compFaceIntegrals
        (
            points,
            faces[cFaces[i]],
            n,
            fW[i],
            A,
            B,
            C,

            Fa,
            Fb,
            Fc,
            Faa,
            Fbb,
            Fcc,
            Faaa,
            Fbbb,
            Fccc,
            Faab,
            Fbbc,
            Fcca
        );

        T0 += n[0] * ((A == 0) ? Fa : ((B == 0) ? Fb : Fc));

        T1[A] += n[A] * Faa;
        T1[B] += n[B] * Fbb;
        T1[C] += n[C] * Fcc;

        T2[A] += n[A] * Faaa;
        T2[B] += n[B] * Fbbb;
        T2[C] += n[C] * Fccc;

        TP[A] += n[A] * Faab;
        TP[B] += n[B] * Fbbc;
        TP[C] += n[C] * Fcca;
    }

    T1 /= 2;
    T2 /= 3;
    TP /= 2;
}


// Calculate
// - r: centre of mass
// - J: inertia around origin (point 0,0,0)
void Foam::momentOfIntertia
(
    const pointField& points,
    const faceList& faces,
    const cell& cFaces,
    point& r,
    tensor& J
)
{
    // Face normals
    vectorField fNorm(cFaces.size());
    scalarField fW(cFaces.size());

    forAll(cFaces, i)
    {
        label faceI = cFaces[i];

        const face& f = faces[faceI];

        fNorm[i] = f.normal(points);
        fNorm[i] /= mag(fNorm[i]) + VSMALL;

        fW[i] = - (fNorm[i] & points[f[0]]);
    }


    scalar T0;
    vector T1, T2, TP;

    compVolumeIntegrals
    (
        points,
        faces,
        cFaces,
        fNorm,
        fW,

        T0,
        T1,
        T2,
        TP
    );

    const scalar density = 1.0;  /* assume unit density */

    scalar mass = density * T0;

    /* compute center of mass */
    r = T1 / T0;

    /* compute inertia tensor */
    J.xx() = density * (T2[1] + T2[2]);
    J.yy() = density * (T2[2] + T2[0]);
    J.zz() = density * (T2[0] + T2[1]);
    J.xy() = J.yx() = - density * TP[0];
    J.yz() = J.zy() = - density * TP[1];
    J.zx() = J.xz() = - density * TP[2];

    ///* translate inertia tensor to center of mass */
    //J[XX] -= mass * (r[1]*r[1] + r[2]*r[2]);
    //J[YY] -= mass * (r[2]*r[2] + r[0]*r[0]);
    //J[ZZ] -= mass * (r[0]*r[0] + r[1]*r[1]);
    //J[XY] = J[YX] += mass * r[0] * r[1];
    //J[YZ] = J[ZY] += mass * r[1] * r[2];
    //J[ZX] = J[XZ] += mass * r[2] * r[0];
}



// ************************************************************************* //
