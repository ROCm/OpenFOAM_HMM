/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "symmTensor.H"
#include "cubicEqn.H"
#include "mathematicalConstants.H"

using namespace Foam::constant::mathematical;


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const Foam::symmTensor::vsType::typeName = "symmTensor";

template<>
const char* const Foam::symmTensor::vsType::componentNames[] =
{
    "xx", "xy", "xz",
          "yy", "yz",
                "zz"
};

template<>
const Foam::symmTensor Foam::symmTensor::vsType::vsType::zero
(
    symmTensor::uniform(0)
);

template<>
const Foam::symmTensor Foam::symmTensor::vsType::one
(
    symmTensor::uniform(1)
);

template<>
const Foam::symmTensor Foam::symmTensor::vsType::max
(
    symmTensor::uniform(VGREAT)
);

template<>
const Foam::symmTensor Foam::symmTensor::vsType::min
(
    symmTensor::uniform(-VGREAT)
);

template<>
const Foam::symmTensor Foam::symmTensor::vsType::rootMax
(
    symmTensor::uniform(ROOTVGREAT)
);

template<>
const Foam::symmTensor Foam::symmTensor::vsType::rootMin
(
    symmTensor::uniform(-ROOTVGREAT)
);

template<>
const Foam::symmTensor Foam::symmTensor::I
(
    1, 0, 0,
       1, 0,
          1
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::vector Foam::eigenValues(const symmTensor& t)
{
    // Coefficients of the characteristic cubic polynomial (a = 1)
    const scalar b =
      - t.xx() - t.yy() - t.zz();
    const scalar c =
        t.xx()*t.yy() + t.xx()*t.zz() + t.yy()*t.zz()
      - t.xy()*t.yx() - t.yz()*t.zy() - t.zx()*t.xz();
    const scalar d =
      - t.xx()*t.yy()*t.zz()
      - t.xy()*t.yz()*t.zx() - t.xz()*t.zy()*t.yx()
      + t.xx()*t.yz()*t.zy() + t.yy()*t.zx()*t.xz() + t.zz()*t.xy()*t.yx();

    // Determine the roots of the characteristic cubic polynomial
    Roots<3> roots(cubicEqn(1, b, c, d).roots());

    vector eVals(Zero);

    // Check the root types
    forAll(roots, i)
    {
        switch (roots.type(i))
        {
            case roots::real:
                eVals[i] = roots[i];
                break;
            case roots::complex:
            case roots::posInf:
            case roots::negInf:
            case roots::nan:
                FatalErrorInFunction
                    << "Eigenvalue computation fails for symmTensor: " << t
                    << "due to the non-real root = " << roots[i]
                    << exit(FatalError);
        }
    }

    // Sort the eigenvalues into ascending order
    if (eVals.x() > eVals.y())
    {
        Swap(eVals.x(), eVals.y());
    }
    if (eVals.y() > eVals.z())
    {
        Swap(eVals.y(), eVals.z());
    }
    if (eVals.x() > eVals.y())
    {
        Swap(eVals.x(), eVals.y());
    }

    return eVals;
}


Foam::vector Foam::eigenVector
(
    const symmTensor& T,
    const scalar eVal,
    const vector& standardBasis1,
    const vector& standardBasis2
)
{
    // Construct the characteristic equation system for this eigenvalue
    const tensor A(T - eVal*I);

    // Determinants of the 2x2 sub-matrices used to find the eigenvectors
    // Sub-determinants for a unique eigenvenvalue
    scalar sd0 = A.yy()*A.zz() - A.yz()*A.zy();
    scalar sd1 = A.zz()*A.xx() - A.zx()*A.xz();
    scalar sd2 = A.xx()*A.yy() - A.xy()*A.yx();
    scalar magSd0 = mag(sd0);
    scalar magSd1 = mag(sd1);
    scalar magSd2 = mag(sd2);

    // Evaluate the eigenvector using the largest sub-determinant
    if (magSd0 >= magSd1 && magSd0 >= magSd2 && magSd0 > SMALL)
    {
        const vector eVec
        (
            1,
            (A.yz()*A.zx() - A.zz()*A.yx())/sd0,
            (A.zy()*A.yx() - A.yy()*A.zx())/sd0
        );

        #ifdef FULLDEBUG
        if (mag(eVec) < SMALL)
        {
            FatalErrorInFunction
                << "Eigenvector magnitude should be non-zero:"
                << "mag(eigenvector) = " << mag(eVec)
                << abort(FatalError);
        }
        #endif

        return eVec/mag(eVec);
    }
    else if (magSd1 >= magSd2 && magSd1 > SMALL)
    {
        const vector eVec
        (
            (A.xz()*A.zy() - A.zz()*A.xy())/sd1,
            1,
            (A.zx()*A.xy() - A.xx()*A.zy())/sd1
        );

        #ifdef FULLDEBUG
        if (mag(eVec) < SMALL)
        {
            FatalErrorInFunction
                << "Eigenvector magnitude should be non-zero:"
                << "mag(eigenvector) = " << mag(eVec)
                << abort(FatalError);
        }
        #endif

        return eVec/mag(eVec);
    }
    else if (magSd2 > SMALL)
    {
        const vector eVec
        (
            (A.xy()*A.yz() - A.yy()*A.xz())/sd2,
            (A.yx()*A.xz() - A.xx()*A.yz())/sd2,
            1
        );

        #ifdef FULLDEBUG
        if (mag(eVec) < SMALL)
        {
            FatalErrorInFunction
                << "Eigenvector magnitude should be non-zero:"
                << "mag(eigenvector) = " << mag(eVec)
                << abort(FatalError);
        }
        #endif

        return eVec/mag(eVec);
    }

    // Sub-determinants for a repeated eigenvalue
    sd0 = A.yy()*standardBasis1.z() - A.yz()*standardBasis1.y();
    sd1 = A.zz()*standardBasis1.x() - A.zx()*standardBasis1.z();
    sd2 = A.xx()*standardBasis1.y() - A.xy()*standardBasis1.x();
    magSd0 = mag(sd0);
    magSd1 = mag(sd1);
    magSd2 = mag(sd2);

    // Evaluate the eigenvector using the largest sub-determinant
    if (magSd0 >= magSd1 && magSd0 >= magSd2 && magSd0 > SMALL)
    {
        const vector eVec
        (
            1,
            (A.yz()*standardBasis1.x() - standardBasis1.z()*A.yx())/sd0,
            (standardBasis1.y()*A.yx() - A.yy()*standardBasis1.x())/sd0
        );

        #ifdef FULLDEBUG
        if (mag(eVec) < SMALL)
        {
            FatalErrorInFunction
                << "Eigenvector magnitude should be non-zero:"
                << "mag(eigenvector) = " << mag(eVec)
                << abort(FatalError);
        }
        #endif

        return eVec/mag(eVec);
    }
    else if (magSd1 >= magSd2 && magSd1 > SMALL)
    {
        const vector eVec
        (
            (standardBasis1.z()*A.zy() - A.zz()*standardBasis1.y())/sd1,
            1,
            (A.zx()*standardBasis1.y() - standardBasis1.x()*A.zy())/sd1
        );

        #ifdef FULLDEBUG
        if (mag(eVec) < SMALL)
        {
            FatalErrorInFunction
                << "Eigenvector magnitude should be non-zero:"
                << "mag(eigenvector) = " << mag(eVec)
                << abort(FatalError);
        }
        #endif

        return eVec/mag(eVec);
    }
    else if (magSd2 > SMALL)
    {
        const vector eVec
        (
            (A.xy()*standardBasis1.z() - standardBasis1.y()*A.xz())/sd2,
            (standardBasis1.x()*A.xz() - A.xx()*standardBasis1.z())/sd2,
            1
        );

        #ifdef FULLDEBUG
        if (mag(eVec) < SMALL)
        {
            FatalErrorInFunction
                << "Eigenvector magnitude should be non-zero:"
                << "mag(eigenvector) = " << mag(eVec)
                << abort(FatalError);
        }
        #endif

        return eVec/mag(eVec);
    }

    // Triple eigenvalue
    return standardBasis1^standardBasis2;
}


Foam::tensor Foam::eigenVectors
(
    const symmTensor& T,
    const vector& eVals
)
{
    vector Ux(1, 0, 0), Uy(0, 1, 0), Uz(0, 0, 1);

    Ux = eigenVector(T, eVals.x(), Uy, Uz);
    Uy = eigenVector(T, eVals.y(), Uz, Ux);
    Uz = eigenVector(T, eVals.z(), Ux, Uy);

    return tensor(Ux, Uy, Uz);
}


Foam::tensor Foam::eigenVectors(const symmTensor& T)
{
    const vector eVals(eigenValues(T));

    return eigenVectors(T, eVals);
}


// ************************************************************************* //
