/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2017 OpenFOAM Foundation
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
    Foam::cylinderToFace

Description
    A topoSetFaceSource to select faces based on face centres inside a cylinder.

    \heading Dictionary parameters
    \table
        Property    | Description                       | Required | Default
        p1          | Coordinate of endpoint            | yes   |
        p2          | Coordinate of endpoint            | yes   |
        radius      | Cylinder (outer) radius           | yes   |
        innerRadius | Cylinder inner radius             | no    | 0
    \endtable

SourceFiles
    cylinderToFace.C

\*---------------------------------------------------------------------------*/

#ifndef cylinderToFace_H
#define cylinderToFace_H

#include "topoSetFaceSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class cylinderToFace Declaration
\*---------------------------------------------------------------------------*/

class cylinderToFace
:
    public topoSetFaceSource
{

    // Private data

        //- Add usage string
        static addToUsageTable usage_;

        //- First point on cylinder axis
        point point1_;

        //- Second point on cylinder axis
        point point2_;

        //- Outer radius
        scalar radius_;

        //- Inner radius
        scalar innerRadius_;


    // Private Member Functions

        void combine(topoSet& set, const bool add) const;


public:

    //- Runtime type information
    TypeName("cylinderToFace");


    // Constructors

        //- Construct from components
        cylinderToFace
        (
            const polyMesh& mesh,
            const point& point1,
            const point& point2,
            const scalar radius,
            const scalar innerRadius = 0
        );

        //- Construct from dictionary
        cylinderToFace(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        cylinderToFace(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~cylinderToFace() = default;


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
