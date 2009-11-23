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

// calcDistances generates the distances_ lookup table (cumulative
// distance along the line) from the individual vectors to the points

void Foam::polyLine::calcDistances()
{
    distances_.setSize(controlPoints_.size());

    if (distances_.size())
    {
        distances_[0] = 0.0;

        for (label i=1; i < distances_.size(); i++)
        {
            distances_[i] = distances_[i-1] +
                mag(controlPoints_[i] - controlPoints_[i-1]);
        }

        // normalize on the interval 0-1
        lineLength_ = distances_[distances_.size()-1];
        for (label i=1; i < distances_.size() - 1; i++)
        {
            distances_[i] /= lineLength_;
        }
        distances_[distances_.size()-1] = 1.0;
    }
    else
    {
        lineLength_ = 0.0;
    }
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polyLine::polyLine(const pointField& ps)
:
    controlPoints_(ps),
    distances_(0),
    lineLength_(0.0)
{
    calcDistances();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::pointField& Foam::polyLine::controlPoints() const
{
    return controlPoints_;
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
        return controlPoints_[0];
    }
    else if (lambda > 1 - SMALL)
    {
        return controlPoints_[controlPoints_.size()-1];
    }


    // search table of cumulative distances to find which line-segment
    // we are on. Check the upper bound.

    label i = 1;
    while (distances_[i] < lambda)
    {
        i++;
    }
    i--;   // we now want the lower bound


    // linear interpolation
    return
    (
        controlPoints_[i]
      + ( controlPoints_[i+1] - controlPoints_[i] )
      * ( lambda - distances_[i] ) / ( distances_[i+1] - distances_[i] )
    );
}


Foam::scalar Foam::polyLine::length() const
{
    return lineLength_;
}


// ************************************************************************* //
