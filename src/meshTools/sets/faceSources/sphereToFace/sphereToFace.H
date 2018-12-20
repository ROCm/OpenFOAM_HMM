/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::sphereToFace

Description
    A topoSetFaceSource to select faces based on faces centres inside sphere.

    \heading Dictionary parameters
    \table
        Property    | Description                           | Required | Default
        origin      | The origin (centre) of the sphere     | yes   |
        radius      | The (outer) radius of sphere          | yes   |
        innerRadius | The inner radius of sphere            | no    | 0
        centre      | Alternative for 'origin'              | no    |
    \endtable

SourceFiles
    sphereToFace.C

\*---------------------------------------------------------------------------*/

#ifndef sphereToFace_H
#define sphereToFace_H

#include "topoSetFaceSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class sphereToFace Declaration
\*---------------------------------------------------------------------------*/

class sphereToFace
:
    public topoSetFaceSource
{

    // Private data

        //- Add usage string
        static addToUsageTable usage_;

        //- Centre point of the sphere
        point origin_;

        //- The outer radius of the sphere
        scalar radius_;

        //- The inner radius of the sphere
        scalar innerRadius_;


    // Private Member Functions

        void combine(topoSet& set, const bool add) const;


public:

    //- Runtime type information
    TypeName("sphereToFace");


    // Constructors

        //- Construct from components
        sphereToFace
        (
            const polyMesh& mesh,
            const point& origin,
            const scalar radius,
            const scalar innerRadius = 0
        );

        //- Construct from dictionary
        sphereToFace(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        sphereToFace(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~sphereToFace() = default;


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
