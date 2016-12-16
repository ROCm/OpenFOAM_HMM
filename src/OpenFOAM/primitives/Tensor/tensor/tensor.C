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

\*---------------------------------------------------------------------------*/

#include "tensor.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const Foam::tensor::vsType::typeName = "tensor";

template<>
const char* const Foam::tensor::vsType::componentNames[] =
{
    "xx", "xy", "xz",
    "yx", "yy", "yz",
    "zx", "zy", "zz"
};

template<>
const Foam::tensor Foam::tensor::vsType::zero(tensor::uniform(0));

template<>
const Foam::tensor Foam::tensor::vsType::one(tensor::uniform(1));

template<>
const Foam::tensor Foam::tensor::vsType::max(tensor::uniform(VGREAT));

template<>
const Foam::tensor Foam::tensor::vsType::min(tensor::uniform(-VGREAT));

template<>
const Foam::tensor Foam::tensor::vsType::rootMax(tensor::uniform(ROOTVGREAT));

template<>
const Foam::tensor Foam::tensor::vsType::rootMin(tensor::uniform(-ROOTVGREAT));

template<>
const Foam::tensor Foam::tensor::I
(
    1, 0, 0,
    0, 1, 0,
    0, 0, 1
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::vector Foam::eigenValues(const tensor& T)
{
    // The eigenvalues
    scalar i, ii, iii;

    // diagonal matrix
    const scalar onDiagMagSum =
        (
            mag(T.xx()) + mag(T.yy()) + mag(T.zz())
        );

    const scalar offDiagMagSum =
        (
            mag(T.xy()) + mag(T.xz()) + mag(T.yx())
          + mag(T.yz()) + mag(T.zx()) + mag(T.zy())
        );

    const scalar magSum = onDiagMagSum + offDiagMagSum;

    if (offDiagMagSum < max(VSMALL, SMALL*magSum))
    {
        i = T.xx();
        ii = T.yy();
        iii = T.zz();
    }

    // non-diagonal matrix
    else
    {
        // Coefficients of the characteristic polynmial
        // x^3 + a*x^2 + b*x + c = 0
        scalar a =
           - T.xx() - T.yy() - T.zz();

        scalar b =
            T.xx()*T.yy() + T.xx()*T.zz() + T.yy()*T.zz()
          - T.xy()*T.yx() - T.yz()*T.zy() - T.zx()*T.xz();

        scalar c =
          - T.xx()*T.yy()*T.zz()
          - T.xy()*T.yz()*T.zx() - T.xz()*T.zy()*T.yx()
          + T.xx()*T.yz()*T.zy() + T.yy()*T.zx()*T.xz() + T.zz()*T.xy()*T.yx();

        // Auxillary variables
        scalar aBy3 = a/3;

        scalar P = (a*a - 3*b)/9; // == -p_wikipedia/3
        scalar PPP = P*P*P;

        scalar Q = (2*a*a*a - 9*a*b + 27*c)/54; // == q_wikipedia/2
        scalar QQ = Q*Q;

        // Three identical roots
        if (mag(P) < SMALL*sqr(magSum) && mag(Q) < SMALL*pow3(magSum))
        {
            return vector(- aBy3, - aBy3, - aBy3);
        }

        // Two identical roots and one distinct root
        else if (mag(PPP - QQ) < SMALL*pow6(magSum))
        {
            scalar sqrtP = sqrt(P);
            scalar signQ = sign(Q);

            i = ii = signQ*sqrtP - aBy3;
            iii = - 2*signQ*sqrtP - aBy3;
        }

        // Three distinct roots
        else if (PPP > QQ)
        {
            scalar sqrtP = sqrt(P);
            scalar value = cos(acos(Q/sqrt(PPP))/3);
            scalar delta = sqrt(3 - 3*value*value);

            i = - 2*sqrtP*value - aBy3;
            ii = sqrtP*(value + delta) - aBy3;
            iii = sqrtP*(value - delta) - aBy3;
        }

        // One real root, two imaginary roots
        // based on the above logic, PPP must be less than QQ
        else
        {
            WarningInFunction
                << "complex eigenvalues detected for tensor: " << T
                << endl;

            if (mag(P) < SMALL*sqr(magSum))
            {
                i = cbrt(QQ/2);
            }
            else
            {
                scalar w = cbrt(- Q - sqrt(QQ - PPP));
                i = w + P/w - aBy3;
            }

            return vector(-VGREAT, i, VGREAT);
        }
    }

    // Sort the eigenvalues into ascending order
    if (i > ii)
    {
        Swap(i, ii);
    }

    if (ii > iii)
    {
        Swap(ii, iii);
    }

    if (i > ii)
    {
        Swap(i, ii);
    }

    return vector(i, ii, iii);
}


Foam::vector Foam::eigenVector
(
    const tensor& T,
    const scalar lambda,
    const vector& direction1,
    const vector& direction2
)
{
    // Construct the linear system for this eigenvalue
    tensor A(T - lambda*I);

    // Determinants of the 2x2 sub-matrices used to find the eigenvectors
    scalar sd0, sd1, sd2;
    scalar magSd0, magSd1, magSd2;

    // Sub-determinants for a unique eivenvalue
    sd0 = A.yy()*A.zz() - A.yz()*A.zy();
    sd1 = A.zz()*A.xx() - A.zx()*A.xz();
    sd2 = A.xx()*A.yy() - A.xy()*A.yx();
    magSd0 = mag(sd0);
    magSd1 = mag(sd1);
    magSd2 = mag(sd2);

    // Evaluate the eigenvector using the largest sub-determinant
    if (magSd0 >= magSd1 && magSd0 >= magSd2 && magSd0 > SMALL)
    {
        vector ev
        (
            1,
            (A.yz()*A.zx() - A.zz()*A.yx())/sd0,
            (A.zy()*A.yx() - A.yy()*A.zx())/sd0
        );

        return ev/mag(ev);
    }
    else if (magSd1 >= magSd2 && magSd1 > SMALL)
    {
        vector ev
        (
            (A.xz()*A.zy() - A.zz()*A.xy())/sd1,
            1,
            (A.zx()*A.xy() - A.xx()*A.zy())/sd1
        );

        return ev/mag(ev);
    }
    else if (magSd2 > SMALL)
    {
        vector ev
        (
            (A.xy()*A.yz() - A.yy()*A.xz())/sd2,
            (A.yx()*A.xz() - A.xx()*A.yz())/sd2,
            1
        );

        return ev/mag(ev);
    }

    // Sub-determinants for a repeated eigenvalue
    sd0 = A.yy()*direction1.z() - A.yz()*direction1.y();
    sd1 = A.zz()*direction1.x() - A.zx()*direction1.z();
    sd2 = A.xx()*direction1.y() - A.xy()*direction1.x();
    magSd0 = mag(sd0);
    magSd1 = mag(sd1);
    magSd2 = mag(sd2);

    // Evaluate the eigenvector using the largest sub-determinant
    if (magSd0 >= magSd1 && magSd0 >= magSd2 && magSd0 > SMALL)
    {
        vector ev
        (
            1,
            (A.yz()*direction1.x() - direction1.z()*A.yx())/sd0,
            (direction1.y()*A.yx() - A.yy()*direction1.x())/sd0
        );

        return ev/mag(ev);
    }
    else if (magSd1 >= magSd2 && magSd1 > SMALL)
    {
        vector ev
        (
            (direction1.z()*A.zy() - A.zz()*direction1.y())/sd1,
            1,
            (A.zx()*direction1.y() - direction1.x()*A.zy())/sd1
        );

        return ev/mag(ev);
    }
    else if (magSd2 > SMALL)
    {
        vector ev
        (
            (A.xy()*direction1.z() - direction1.y()*A.xz())/sd2,
            (direction1.x()*A.xz() - A.xx()*direction1.z())/sd2,
            1
        );

        return ev/mag(ev);
    }

    // Triple eigenvalue
    return direction1^direction2;
}


Foam::tensor Foam::eigenVectors(const tensor& T, const vector& lambdas)
{
    vector Ux(1, 0, 0), Uy(0, 1, 0), Uz(0, 0, 1);

    Ux = eigenVector(T, lambdas.x(), Uy, Uz);
    Uy = eigenVector(T, lambdas.y(), Uz, Ux);
    Uz = eigenVector(T, lambdas.z(), Ux, Uy);

    return tensor(Ux, Uy, Uz);
}


Foam::tensor Foam::eigenVectors(const tensor& T)
{
    const vector lambdas(eigenValues(T));

    return eigenVectors(T, lambdas);
}


Foam::vector Foam::eigenValues(const symmTensor& T)
{
    return eigenValues(tensor(T));
}


Foam::vector Foam::eigenVector
(
    const symmTensor& T,
    const scalar lambda,
    const vector& direction1,
    const vector& direction2
)
{
    return eigenVector(tensor(T), lambda, direction1, direction2);
}


Foam::tensor Foam::eigenVectors(const symmTensor& T, const vector& lambdas)
{
    return eigenVectors(tensor(T), lambdas);
}


Foam::tensor Foam::eigenVectors(const symmTensor& T)
{
    return eigenVectors(tensor(T));
}


// ************************************************************************* //
