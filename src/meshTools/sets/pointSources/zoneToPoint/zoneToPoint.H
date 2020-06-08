/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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
    A \c topoSetPointSource to convert \c pointZone(s) to a \c pointSet.

    Operands:
    \table
      Operand   | Type         | Location
      input     | pointZone(s) | $FOAM_CASE/constant/polyMesh/pointZones
      output    | pointSet     | $FOAM_CASE/constant/polyMesh/sets/\<set\>
    \endtable

Usage
    Minimal example by using \c system/topoSetDict.actions:
    \verbatim
    {
        // Mandatory (inherited) entries
        name        <name>;
        type        pointSet;
        action      <action>;

        // Mandatory entries
        source      zoneToPoint;

        // Conditional mandatory entries
        // Select either of the below

        // Option-1
        zones
        (
            <pointZoneName0>
            <pointZoneName1>
            ...
        );

        // Option-2
        zone     <pointZoneName>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                           | Type | Req'd | Dflt
      name       | Name of pointSet                      | word |  yes  | -
      type       | Type name: pointSet                   | word |  yes  | -
      action     | Action applied on points - see below  | word |  yes  | -
      source     | Source name: zoneToPoint              | word |  yes  | -
    \endtable

    Options for the \c action entry:
    \verbatim
      new      | Create a new pointSet from selected points
      add      | Add selected points into this pointSet
      subtract | Remove selected points from this pointSet
    \endverbatim

    Options for the conditional mandatory entries:
    \verbatim
      Entry    | Description                | Type     | Req'd  | Dflt
      zones    | Names of input pointZones  | wordList | cond'l | -
      zone     | Name of input pointZone    | word     | cond'l | -
    \verbatim

Note
    The order of precedence among the conditional mandatory entries from the
    highest to the lowest is \c zones, and \c zone.

See also
    - Foam::topoSetSource
    - Foam::topoSetPointSource

SourceFiles
    zoneToPoint.C

\*---------------------------------------------------------------------------*/

#ifndef zoneToPoint_H
#define zoneToPoint_H

#include "topoSetPointSource.H"
#include "wordRes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class zoneToPoint Declaration
\*---------------------------------------------------------------------------*/

class zoneToPoint
:
    public topoSetPointSource
{
    // Private Data

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
