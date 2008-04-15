/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
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
    Foam::boundBox

Description
    A bounding box defined in terms of the points at its extremities.

\*---------------------------------------------------------------------------*/

#ifndef boundBox_H
#define boundBox_H

#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class boundBox;

Ostream& operator<<(Ostream& os, const boundBox& b);


/*---------------------------------------------------------------------------*\
                           Class boundBox Declaration
\*---------------------------------------------------------------------------*/

class boundBox
{
    // Private data

        //- Minimum and maximum describing the bounding box
        point min_, max_;


public:

    // Constructors

        //- Construct null setting points to zero
        boundBox()
        :
            min_(vector::zero),
            max_(vector::zero)
        {}

        //- Construct from components
        boundBox(const point& min, const point& max)
        :
            min_(min),
            max_(max)
        {}

        //- Construct as the bounding box of the given pointField. Does
        //  parallel communication (doReduce = true)
        boundBox(const pointField& points, const bool doReduce = true);

        //- Construct from Istream
        boundBox(Istream&);


    // Member functions

        // Access

            const point& min() const
            {
                return min_;
            }

            const point& max() const
            {
                return max_;
            }

            point& min()
            {
                return min_;
            }

            point& max()
            {
                return max_;
            }

            scalar mag() const
            {
                return ::Foam::mag(max_ - min_);
            }


        // Query

            //- Intersects other boundingbox?
            bool intersects(const boundBox& bb) const
            {
                if
                (
                    (min_.x() <= bb.max().x()) &&
                    (min_.y() <= bb.max().y()) &&
                    (min_.z() <= bb.max().z()) &&

                    (max_.x() >= bb.min().x()) &&
                    (max_.y() >= bb.min().y()) &&
                    (max_.z() >= bb.min().z())
                )
                {
                    return true;
                }
                else
                {
                    return false;
                }
            }


            //- Contains a point?
            bool contains(const point& pt) const
            {
                return
                    pt.x() >= min().x()
                 && pt.y() >= min().y()
                 && pt.z() >= min().z()
                 && pt.x() <= max().x()
                 && pt.y() <= max().y()
                 && pt.z() <= max().z();
            }


    // Ostream operator

        friend Ostream& operator<<(Ostream& os, const boundBox& b);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
