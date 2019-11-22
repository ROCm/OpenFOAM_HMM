/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2018-2019 OpenCFD Ltd.
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
    Foam::faceZoneToFaceZone

Description
    A topoSetSource to select faces based on usage in another faceSet.

    \heading Dictionary parameters
    \table
        Property    | Description                        | Required  | Default
        zone        | Name of input face zone (no regex) | yes   |
    \endtable

SourceFiles
    faceZoneToFaceZone.C

\*---------------------------------------------------------------------------*/

#ifndef faceZoneToFaceZone_H
#define faceZoneToFaceZone_H

#include "topoSetSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class faceZoneToFaceZone Declaration
\*---------------------------------------------------------------------------*/

class faceZoneToFaceZone
:
    public topoSetSource
{
    // Private data

        //- Add usage string
        static addToUsageTable usage_;

        //- Name of set to use
        word setName_;

public:

    //- Runtime type information
    TypeName("faceZoneToFaceZone");

    // Constructors

        //- Construct from components
        faceZoneToFaceZone
        (
            const polyMesh& mesh,
            const word& setName
        );

        //- Construct from dictionary
        faceZoneToFaceZone(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        faceZoneToFaceZone(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~faceZoneToFaceZone() = default;


    // Member Functions

        //- The source category is a faceZone
        virtual topoSetSource::sourceType setType() const
        {
            return FACEZONE_SOURCE;
        }

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
