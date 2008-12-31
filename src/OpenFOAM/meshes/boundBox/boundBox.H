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

    // Static data members

        //- A very large boundBox: min/max == -/+ VGREAT
        static const boundBox greatBox;

        //- A very large inverted boundBox: min/max == +/- VGREAT
        static const boundBox invertedBox;


    // Constructors

        //- Construct null, setting points to zero
        boundBox()
        :
            min_(point::zero),
            max_(point::zero)
        {}

        //- Construct from components
        boundBox(const point& min, const point& max)
        :
            min_(min),
            max_(max)
        {}

        //- Construct as the bounding box of the given pointField.
        //  Does parallel communication (doReduce = true)
        boundBox(const pointField& points, const bool doReduce = true);

        //- Construct from Istream
        boundBox(Istream&);


    // Member functions

        // Access

            //- Minimum describing the bounding box
            const point& min() const
            {
                return min_;
            }

            //- Maximum describing the bounding box
            const point& max() const
            {
                return max_;
            }

            //- Minimum describing the bounding box, non-const access
            point& min()
            {
                return min_;
            }

            //- Maximum describing the bounding box, non-const access
            point& max()
            {
                return max_;
            }

            //- The bounding box span (from minimum to maximum)
            vector span() const
            {
                return (max_ - min_);
            }

            //- The magnitude of the bounding box span
            scalar mag() const
            {
                return ::Foam::mag(max_ - min_);
            }


        // Query

            //- Intersects other boundingBox?
            bool overlaps(const boundBox& bb) const
            {
                return
                (
                    min_.x() <= bb.max().x() && max_.x() >= bb.min().x()
                 && min_.y() <= bb.max().y() && max_.y() >= bb.min().y()
                 && min_.z() <= bb.max().z() && max_.z() >= bb.min().z()
                );
            }

            //- Contains a point?
            bool contains(const point& pt) const
            {
                return
                (
                    pt.x() >= min().x() && pt.x() <= max().x()
                 && pt.y() >= min().y() && pt.y() <= max().y()
                 && pt.z() >= min().z() && pt.z() <= max().z()
                );
            }


    // Friend Operators

        friend bool operator==(const boundBox& a, const boundBox& b)
        {
            return (a.min_ == b.min_) && (a.max_ == b.max_);
        }

        friend bool operator!=(const boundBox& a, const boundBox& b)
        {
            return !(a == b);
        }


    // IOstream operator

        friend Istream& operator>>(Istream& is, boundBox&);
        friend Ostream& operator<<(Ostream& os, const boundBox&);
};


//- Specify data associated with boundBox type are contiguous
template<>
inline bool contiguous<boundBox>() {return contiguous<point>();}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
