/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2011 OpenCFD Ltd.
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

Class
    Foam::densityWeightedStochastic

Description
    Choose random points inside the domain and place them with a probability
    proportional to the target density of points.

SourceFiles
    densityWeightedStochastic.C

\*---------------------------------------------------------------------------*/

#ifndef densityWeightedStochastic_H
#define densityWeightedStochastic_H

#include "initialPointsMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class densityWeightedStochastic Declaration
\*---------------------------------------------------------------------------*/

class densityWeightedStochastic
:
    public initialPointsMethod
{

private:

    // Private data

        //- The total volume to be filled
        scalar totalVolume_;

        //- Working variable for minimum cell size, a starting value may be
        //  specified in the dictionary
        mutable scalar minCellSize_;

        //- Smallest minimum cell size allowed, i.e. to avoid high initial
        //  population of areas of small size
        scalar minCellSizeLimit_;


public:

    //- Runtime type information
    TypeName("densityWeightedStochastic");

    // Constructors

        //- Construct from components
        densityWeightedStochastic
        (
            const dictionary& initialPointsDict,
            const conformalVoronoiMesh& cvMesh
        );


    //- Destructor
    virtual ~densityWeightedStochastic()
    {}


    // Member Functions

        //- Return the initial points for the conformalVoronoiMesh
        virtual std::vector<Vb::Point> initialPoints() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
