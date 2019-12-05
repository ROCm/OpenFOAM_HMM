/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
    Copyright (C) 2017-2018 OpenCFD Ltd.
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

Class
    Foam::coordSet

Description
    Holds list of sampling positions

SourceFiles
    coordSet.C

\*---------------------------------------------------------------------------*/

#ifndef coordSet_H
#define coordSet_H

#include "pointField.H"
#include "word.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class coordSet Declaration
\*---------------------------------------------------------------------------*/

class coordSet
:
    public pointField
{
public:

    // Public data types

        //- Enumeration defining the output format for coordinates
        enum class coordFormat
        {
            XYZ,
            X,
            Y,
            Z,
            DISTANCE
        };

        //- String representation of coordFormat enum
        static const Enum<coordFormat> coordFormatNames;

private:

        //- Check for consistent dimensions of points and curve distance
        void checkDimensions() const;


protected:

        //- Name
        const word name_;

        //- Axis write type
        const coordFormat axis_;

        //- Cumulative distance "distance" write specifier.
        scalarList curveDist_;


public:

    // Constructors

        //- Construct from components
        //  Note: curveDist will be empty
        coordSet(const word& name, const coordFormat axisType);

        //- Construct from components
        //  Note: curveDist will be empty
        coordSet(const word& name, const word& axis);

        //- Copy construct from components
        coordSet
        (
            const word& name,
            const word& axis,
            const List<point>& points,
            const scalarList& curveDist
        );

        //- Move construct from components
        coordSet
        (
            const word& name,
            const word& axis,
            List<point>&& points,
            scalarList&& curveDist
        );


    // Member functions

        const word& name() const
        {
            return name_;
        }

        word axis() const
        {
            return coordFormatNames[axis_];
        }

        //- Set the points
        void setPoints(const List<point>& newPoints)
        {
            static_cast<pointField&>(*this) = newPoints;
        }

        //- Set the points
        void setPoints(List<point>&& newPoints)
        {
            static_cast<pointField&>(*this) = std::move(newPoints);
        }

        //- Return the cumulative distance
        const scalarList& curveDist() const
        {
            return curveDist_;
        }

        //- Set the cumulative distance
        void setCurveDist(const scalarList& curveDist)
        {
            curveDist_ = curveDist;
            checkDimensions();
        }

        //- Set the cumulative distance
        void setCurveDist(scalarList&& curveDist)
        {
            curveDist_ = std::move(curveDist);
            checkDimensions();
        }

        //- Is axis specification a vector
        bool hasVectorAxis() const;

        //- Get coordinate of point according to axis specification.
        //  If axis="distance" is the curveDist[index]
        scalar scalarCoord(const label index) const;

        //- Get point according to axis="xyz" specification
        vector vectorCoord(const label index) const;

        //- Write to stream
        Ostream& write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
