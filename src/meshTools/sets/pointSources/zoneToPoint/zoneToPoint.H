/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::zoneToPoint

Description
    A topoSetSource to select points based on one or more pointZones.

    \heading Dictionary parameters
    \table
        Property    | Description                           | Required | Default
        zone        | The point zone name or regex          | possibly |
        zones       | The point zone names or regexs        | possibly |
        name        | Older specification for 'zone'        | no    |
    \endtable

Note
    Selection of multiple zones has precedence.

SourceFiles
    zoneToPoint.C

\*---------------------------------------------------------------------------*/

#ifndef zoneToPoint_H
#define zoneToPoint_H

#include "topoSetSource.H"
#include "wordRes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class zoneToPoint Declaration
\*---------------------------------------------------------------------------*/

class zoneToPoint
:
    public topoSetSource
{

    // Private data

        //- Add usage string
        static addToUsageTable usage_;

        //- Matcher for zones
        wordRes selectedZones_;


    // Private Member Functions

        void combine(topoSet& set, const bool add) const;


public:

    //- Runtime type information
    TypeName("zoneToPoint");

    // Constructors

        //- Construct from components
        zoneToPoint(const polyMesh& mesh, const wordRe& zoneName);

        //- Construct from dictionary
        zoneToPoint(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        zoneToPoint(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~zoneToPoint() = default;


    // Member Functions

        virtual topoSetSource::sourceType setType() const
        {
            return POINTSETSOURCE;
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
