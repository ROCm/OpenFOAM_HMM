/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::zoneToFace

Description
    A topoSetFaceSource to select faces based on one of more faceZones.

    \heading Dictionary parameters
    \table
        Property    | Description                           | Required | Default
        zones       | The face zone names or regexs         | possibly |
        zone        | The face zone name or regex           | possibly |
        name        | Older specification for 'zone'        | no    |
    \endtable

Note
    Must specify "zones", "zone" or "name" (highest to lowest precedence).

SourceFiles
    zoneToFace.C

\*---------------------------------------------------------------------------*/

#ifndef zoneToFace_H
#define zoneToFace_H

#include "topoSetFaceSource.H"
#include "wordRes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class zoneToFace Declaration
\*---------------------------------------------------------------------------*/

class zoneToFace
:
    public topoSetFaceSource
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
    TypeName("zoneToFace");

    // Constructors

        //- Construct from components
        zoneToFace(const polyMesh& mesh, const wordRe& zoneName);

        //- Construct from dictionary
        zoneToFace(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        zoneToFace(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~zoneToFace() = default;


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
