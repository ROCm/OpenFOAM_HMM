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
    Foam::zoneToCell

Description
    A \c topoSetCellSource to convert \c cellZone(s) to a \c cellSet.

    Operands:
    \table
      Operand   | Type        | Location
      input     | cellZone(s) | $FOAM_CASE/constant/polyMesh/cellZones
      output    | cellSet     | $FOAM_CASE/constant/polyMesh/sets/\<set\>
    \endtable

Usage
    Minimal example by using \c system/topoSetDict.actions:
    \verbatim
    {
        // Mandatory (inherited) entries
        name        <name>;
        type        cellSet;
        action      <action>;

        // Mandatory entries
        source      zoneToCell;

        // Conditional mandatory entries
        // Select either of the below

        // Option-1
        zones
        (
            <cellZoneName0>
            <cellZoneName1>
            ...
        );

        // Option-2
        zone     <cellZoneName>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                         | Type | Req'd | Dflt
      name       | Name of cellSet                     | word |  yes  | -
      type       | Type name: cellSet                  | word |  yes  | -
      action     | Action applied on cells - see below | word |  yes  | -
      source     | Source name: zoneToCell             | word |  yes  | -
    \endtable

    Options for the \c action entry:
    \verbatim
      new      | Create a new cellSet from selected cells
      add      | Add selected cells into this cellsSet
      subtract | Remove selected cells from this cellSet
    \endverbatim

    Options for the conditional mandatory entries:
    \verbatim
      Entry    | Description                | Type     | Req'd  | Dflt
      zones    | Names of input cellZones   | wordList | cond'l | -
      zone     | Name of input cellZone     | word     | cond'l | -
    \verbatim

Note
    The order of precedence among the conditional mandatory entries from the
    highest to the lowest is \c zones, and \c zone.

See also
    - Foam::topoSetSource
    - Foam::topoSetCellSource

SourceFiles
    zoneToCell.C

\*---------------------------------------------------------------------------*/

#ifndef zoneToCell_H
#define zoneToCell_H

#include "topoSetCellSource.H"
#include "wordRes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class zoneToCell Declaration
\*---------------------------------------------------------------------------*/

class zoneToCell
:
    public topoSetCellSource
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
    TypeName("zoneToCell");


    // Constructors

        //- Construct from components
        zoneToCell(const polyMesh& mesh, const wordRe& zoneName);

        //- Construct from dictionary
        zoneToCell(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        zoneToCell(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~zoneToCell() = default;


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
