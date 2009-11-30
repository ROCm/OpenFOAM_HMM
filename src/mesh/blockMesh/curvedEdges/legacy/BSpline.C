/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "BSpline.H"
#include "simpleMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::pointField Foam::BSpline::findKnots
(
    const pointField& allknots,
    const vector& fstend,
    const vector& sndend
)
{
    const label NKnots = allknots.size();

    // set up 1/6 and 2/3 which are the matrix elements throughout most
    // of the matrix

    register const scalar oneSixth = 1.0/6.0;
    register const scalar twoThird = 2.0/3.0;

    simpleMatrix<vector> M(NKnots+2, 0, vector::zero);

    // set up the matrix
    M[0][0] = -0.5*scalar(NKnots - 1);
    M[0][2] =  0.5*scalar(NKnots - 1);

    for (register label i = 1; i <= NKnots; i++)
    {
        M[i][i-1] = oneSixth;
        M[i][i] = twoThird;
        M[i][i+1] = oneSixth;
    }

    M[NKnots+1][NKnots-1] = -0.5*scalar(NKnots - 1);
    M[NKnots+1][NKnots+1] =  0.5*scalar(NKnots - 1);

    // set up the vector
    for (register label i = 1; i <= NKnots; i++)
    {
        M.source()[i] = allknots[i-1];
    }

    // set the gradients at the ends:

    if (mag(fstend) < 1e-8)
    {
        // default : forward differences on the end knots
        M.source()[0] = allknots[1] - allknots[0];
        M.source()[0] /= mag(M.source()[0]);
    }
    else
    {
        // use the gradient vector provided
        M.source()[0] = fstend/mag(fstend);
    }

    if (mag(sndend)<1e-8)
    {
        // default : forward differences on the end knots
        M.source()[NKnots+1] = M.source()[NKnots-1] - M.source()[NKnots];
        M.source()[NKnots+1] /= mag(M.source()[NKnots+1]);
    }
    else
    {
        // use the gradient vector provided
        M.source()[NKnots+1] = sndend/mag(sndend);
    }


    // invert the equation to find the control knots
    return M.solve();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::BSpline::BSpline
(
    const pointField& Knots,
    const vector& fstend,
    const vector& sndend
)
:
    spline(findKnots(Knots, fstend, sndend))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::point Foam::BSpline::realPosition(const scalar mu) const
{
    return spline::position(mu);
}


Foam::point Foam::BSpline::position(const scalar mu) const
{
    return spline::position((1.0/(nKnots() - 1))*(1.0 + mu*(nKnots() - 3)));
}


Foam::scalar Foam::BSpline::length() const
{
    notImplemented("BSpline::length() const");
    return 1.0;
}


// ************************************************************************* //
