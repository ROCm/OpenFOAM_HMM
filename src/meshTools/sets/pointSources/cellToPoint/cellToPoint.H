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
    Foam::cellToPoint

Description
    A \c topoSetPointSource to select all the points from given \c cellSet(s).

    Operands:
    \table
      Operand   | Type       | Location
      input     | cellSet(s) | $FOAM_CASE/constant/polyMesh/sets/\<set\>
      output    | pointSet   | $FOAM_CASE/constant/polyMesh/sets/\<set\>
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
        source      cellToPoint;
        option      <option>;

        // Conditional mandatory entries
        // Select either of the below

        // Option-1
        sets
        (
            <pointSetName1>
            <pointSetName2>
            ...
        );

        // Option-2
        set <pointSetName>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                          | Type | Req'd | Dflt
      name       | Name of pointSet                     | word |  yes  | -
      type       | Type name: pointSet                  | word |  yes  | -
      action     | Action applied on points - see below | word |  yes  | -
      source     | Source name: cellToPoint             | word |  yes  | -
      option     | Selection type - see below           | word |  yes  | -
    \endtable

    Options for the \c action entry:
    \verbatim
      new      | Create a new pointSet from selected points
      add      | Add selected points into this pointSet
      subtract | Remove selected points from this pointSet
    \endverbatim

    Options for the \c option entry:
    \verbatim
      all      | Select all points of cells in the cellSet
    \endverbatim

    Options for the conditional mandatory entries:
    \verbatim
      Entry    | Description                    | Type     | Req'd  | Dflt
      sets     | Names of input cellSets        | wordList | cond'l | -
      set      | Name of input cellSet          | word     | cond'l | -
    \endverbatim

Note
    The order of precedence among the conditional mandatory entries from the
    highest to the lowest is \c sets, and \c set.

See also
    - Foam::topoSetSource
    - Foam::topoSetPointSource

SourceFiles
    cellToPoint.C

\*---------------------------------------------------------------------------*/

#ifndef cellToPoint_H
#define cellToPoint_H

#include "topoSetPointSource.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class cellToPoint Declaration
\*---------------------------------------------------------------------------*/

class cellToPoint
:
    public topoSetPointSource
{

public:
        //- Enumeration defining the valid options
        enum cellAction
        {
            ALL
        };


private:

    // Private Data

        //- Add usage string
        static addToUsageTable usage_;

        static const Enum<cellAction> cellActionNames_;

        //- Names of sets to use
        wordList names_;

        //- Option
        cellAction option_;


    // Private Member Functions

        void combine(topoSet& set, const bool add, const word& setName) const;


public:

    //- Runtime type information
    TypeName("cellToPoint");


    // Constructors

        //- Construct from components
        cellToPoint
        (
            const polyMesh& mesh,
            const word& setName,
            const cellAction option
        );

        //- Construct from dictionary
        cellToPoint(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        cellToPoint(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~cellToPoint() = default;


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
