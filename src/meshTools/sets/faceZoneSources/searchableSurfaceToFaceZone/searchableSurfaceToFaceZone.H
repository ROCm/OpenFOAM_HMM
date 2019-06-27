/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2012 OpenFOAM Foundation
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
    Foam::searchableSurfaceToFaceZone

Description
    A topoSetSource to select faces based on intersection (of cell-cell
    vector) with a surface.

    \heading Dictionary parameters
    \table
        Property    | Description                           | Required | Default
        surfaceType | The searchable surface type           | yes   |
        surfaceName | Name for the IOobject                 | no    | mesh-name
        surface     | Same as 'surfaceType'                 | no    |
        name        | Older specification for 'surfaceName' | no    |
    \endtable

SourceFiles
    searchableSurfaceToFaceZone.C

\*---------------------------------------------------------------------------*/

#ifndef searchableSurfaceToFaceZone_H
#define searchableSurfaceToFaceZone_H

#include "topoSetSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class searchableSurface;

/*---------------------------------------------------------------------------*\
                   Class searchableSurfaceToFaceZone Declaration
\*---------------------------------------------------------------------------*/

class searchableSurfaceToFaceZone
:
    public topoSetSource
{
    // Private data

        //- Add usage string
        static addToUsageTable usage_;

        //- Surface
        autoPtr<searchableSurface> surfacePtr_;


protected:

    // Protected Member Functions

        //- Retrieve surface name from dictionary entry
        static word getSurfaceName
        (
            const dictionary& dict,
            const word& defaultName
        );


public:

    //- Runtime type information
    TypeName("searchableSurfaceToFaceZone");

    // Constructors

        //- Construct surface-type from dictionary
        searchableSurfaceToFaceZone
        (
            const word& surfaceType,
            const polyMesh& mesh,
            const dictionary& dict
        );


        //- Construct from dictionary
        searchableSurfaceToFaceZone
        (
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~searchableSurfaceToFaceZone() = default;


    // Member Functions

        virtual topoSetSource::sourceType setType() const
        {
            return FACEZONESOURCE;
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
