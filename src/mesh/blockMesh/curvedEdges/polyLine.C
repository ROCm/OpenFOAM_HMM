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
#include "polyLine.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::polyLine::calcParam()
{
    param_.setSize(points_.size());

    if (param_.size())
    {
        param_[0] = 0.0;

        for (label i=1; i < param_.size(); i++)
        {
            param_[i] = param_[i-1] + mag(points_[i] - points_[i-1]);
        }

        // normalize on the interval 0-1
        lineLength_ = param_.last();
        for (label i=1; i < param_.size() - 1; i++)
        {
            param_[i] /= lineLength_;
        }
        param_.last() = 1.0;
    }
    else
    {
        lineLength_ = 0.0;
    }
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polyLine::polyLine(const pointField& ps)
:
    points_(ps),
    lineLength_(0.0),
    param_(0)
{
    calcParam();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::pointField& Foam::polyLine::points() const
{
    return points_;
}


Foam::label Foam::polyLine::nSegments() const
{
    return points_.size()-1;
}


Foam::label Foam::polyLine::localParameter(scalar& lambda) const
{
    // check range of lambda
    if (lambda < 0 || lambda > 1)
    {
        FatalErrorIn("polyLine::localParameter(scalar&)")
            << "Parameter out-of-range, "
            << "lambda = " << lambda
            << abort(FatalError);
    }

    // check endpoints
    if (lambda < SMALL)
    {
        lambda = 0;
        return 0;
    }
    else if (lambda > 1 - SMALL)
    {
        lambda = 1;
        return nSegments();
    }

    // search table of cumulative distances to find which line-segment
    // we are on. Check the upper bound.

    label segmentI = 1;
    while (param_[segmentI] < lambda)
    {
        segmentI++;
    }
    segmentI--;   // we want the corresponding lower bound

    // the local parameter [0-1] on this line segment
    lambda =
    (
        ( lambda - param_[segmentI] )
      / ( param_[segmentI+1] - param_[segmentI] )
    );

    return segmentI;
}


Foam::point Foam::polyLine::position(const scalar lambda) const
{
    // check range of lambda
    if (lambda < 0 || lambda > 1)
    {
        FatalErrorIn("polyLine::position(const scalar)")
            << "Parameter out of range, "
            << "lambda = " << lambda
            << abort(FatalError);
    }

    // check endpoints
    if (lambda < SMALL)
    {
        return points_[0];
    }
    else if (lambda > 1 - SMALL)
    {
        return points_.last();
    }


    // search table of cumulative distances to find which line-segment
    // we are on. Check the upper bound.

    label segmentI = 1;
    while (param_[segmentI] < lambda)
    {
        ++segmentI;
    }
    --segmentI;   // we now want the lower bound


    // linear interpolation
    return
    (
        points_[segmentI]
      + ( points_[segmentI+1] - points_[segmentI] )
      * ( lambda - param_[segmentI] )
      / ( param_[segmentI+1] - param_[segmentI] )
    );
}


Foam::scalar Foam::polyLine::length() const
{
    return lineLength_;
}


// ************************************************************************* //
