/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::nearestToPoint

Description
    A topoSetPointSource to select points nearest to points.

    \heading Dictionary parameters
    \table
        Property    | Description                       | Required  | Default
        points      | List of selection points          | yes   |
    \endtable

SourceFiles
    nearestToPoint.C

\*---------------------------------------------------------------------------*/

#ifndef nearestToPoint_H
#define nearestToPoint_H

#include "topoSetPointSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class nearestToPoint Declaration
\*---------------------------------------------------------------------------*/

class nearestToPoint
:
    public topoSetPointSource
{

    // Private data

        //- Add usage string
        static addToUsageTable usage_;

        //- Points to select nearest to
        pointField points_;


    // Private Member Functions

        void combine(topoSet& set, const bool add) const;


public:

    //- Runtime type information
    TypeName("nearestToPoint");

    // Constructors

        //- Construct from components, copying points
        nearestToPoint(const polyMesh& mesh, const pointField& points);

        //- Construct from components, moving points
        nearestToPoint(const polyMesh& mesh, pointField&& points);

        //- Construct from dictionary
        nearestToPoint(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        nearestToPoint(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~nearestToPoint() = default;


    // Member Functions

        virtual void applyToSet
        (
            const topoSetSource::setAction action,
            topoSet& set
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
