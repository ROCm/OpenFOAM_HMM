/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2013 OpenFOAM Foundation
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
    Foam::setAndNormalToFaceZone

Description
    A topoSetSource to select faces based on usage in a faceSet, where the
    normal vector is used to orient the faces.

    \heading Dictionary parameters
    \table
        Property    | Description                       | Required  | Default
        faceSet     | Name of input faceSet             | yes   |
        normal      | The normal vector for orientation | yes   |
    \endtable

SourceFiles
    setAndNormalToFaceZone.C

\*---------------------------------------------------------------------------*/

#ifndef setAndNormalToFaceZone_H
#define setAndNormalToFaceZone_H

#include "topoSetSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class setAndNormalToFaceZone Declaration
\*---------------------------------------------------------------------------*/

class setAndNormalToFaceZone
:
    public topoSetSource
{
    // Private data

        //- Add usage string
        static addToUsageTable usage_;

        //- Name of set to use
        word setName_;

        //- Normal used to orient the faces
        vector normal_;


public:

    //- Runtime type information
    TypeName("setAndNormalToFaceZone");

    // Constructors

        //- Construct from components
        setAndNormalToFaceZone
        (
            const polyMesh& mesh,
            const word& setName,
            const vector& normal
        );

        //- Construct from dictionary
        setAndNormalToFaceZone(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        setAndNormalToFaceZone(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~setAndNormalToFaceZone() = default;


    // Member Functions

        virtual sourceType setType() const
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
