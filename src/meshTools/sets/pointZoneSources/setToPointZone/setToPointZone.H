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
    Foam::setToPointZone

Description
    A topoSetSource to select points based on usage in a pointSet.

    \heading Dictionary parameters
    \table
        Property    | Description                       | Required  | Default
        set         | Name of input pointSet            | yes   |
    \endtable

SourceFiles
    setToPointZone.C

\*---------------------------------------------------------------------------*/

#ifndef setToPointZone_H
#define setToPointZone_H

#include "topoSetSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class setToPointZone Declaration
\*---------------------------------------------------------------------------*/

class setToPointZone
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
    TypeName("setToPointZone");

    // Constructors

        //- Construct from components
        setToPointZone
        (
            const polyMesh& mesh,
            const word& setName
        );

        //- Construct from dictionary
        setToPointZone(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        setToPointZone(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~setToPointZone() = default;


    // Member Functions

        //- The source category is a pointZone
        virtual topoSetSource::sourceType setType() const
        {
            return POINTZONE_SOURCE;
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
