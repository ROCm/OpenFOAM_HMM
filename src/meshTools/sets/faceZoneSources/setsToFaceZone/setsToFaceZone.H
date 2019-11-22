/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::setsToFaceZone

Description
    A topoSetSource to select faces based on usage in a faceSet and cellSet

    \heading Dictionary parameters
    \table
        Property    | Description                           | Required | Default
        cellSet     | Name of input cellSet                 | yes   |
        faceSet     | Name of input faceSet                 | yes   |
        flip        | Force flip of faces                   | no    | false
    \endtable

SourceFiles
    setsToFaceZone.C

\*---------------------------------------------------------------------------*/

#ifndef setsToFaceZone_H
#define setsToFaceZone_H

#include "topoSetSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class setsToFaceZone Declaration
\*---------------------------------------------------------------------------*/

class setsToFaceZone
:
    public topoSetSource
{
    // Private data

        //- Add usage string
        static addToUsageTable usage_;

        //- Name of set to use
        const word faceSetName_;

        //- Name of set to use
        const word cellSetName_;

        //- Whether cellSet is slave cells or master cells
        const bool flip_;

public:

    //- Runtime type information
    TypeName("setsToFaceZone");

    // Constructors

        //- Construct from components
        setsToFaceZone
        (
            const polyMesh& mesh,
            const word& faceSetName,
            const word& cellSetName,
            const bool flip
        );

        //- Construct from dictionary
        setsToFaceZone(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        setsToFaceZone(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~setsToFaceZone() = default;


    // Member Functions

        //- The source category is a pointSet
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
