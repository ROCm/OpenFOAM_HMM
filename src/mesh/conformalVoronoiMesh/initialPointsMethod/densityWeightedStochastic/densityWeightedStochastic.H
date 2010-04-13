/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2009-2010 OpenCFD Ltd.
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

SourceFiles
    densityWeightedStochastic.C

\*---------------------------------------------------------------------------*/

#ifndef densityWeightedStochastic_H
#define densityWeightedStochastic_H

#include "Switch.H"
#include "Random.H"
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

        //- 1/(minimum cell size)^3, gives the maximum permissible point
        //  density
        scalar maxDensity_;


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
