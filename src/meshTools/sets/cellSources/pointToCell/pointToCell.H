/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::pointToCell

Description
    A topoSetCellSource to select cells based on usage of points.

    \heading Dictionary parameters
    \table
        Property    | Description                       | Required  | Default
        sets        | The pointSet names to use         | possibly |
        set         | The pointSet name to use          | possibly |
        option      | Selection type (any/edge)         | yes   |
    \endtable

Note
    Must specify "sets" or "set" (highest to lowest precedence).

SourceFiles
    pointToCell.C

\*---------------------------------------------------------------------------*/

#ifndef pointToCell_H
#define pointToCell_H

#include "topoSetCellSource.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class pointToCell Declaration
\*---------------------------------------------------------------------------*/

class pointToCell
:
    public topoSetCellSource
{
public:
        //- Enumeration defining the valid options
        enum pointAction
        {
            ANY,    // Cells using any point in set
            EDGE    // Cells using an edge with both points in set
            //ALL   // Possible extension: cells whose all points are in set
        };

private:

        //- Add usage string
        static addToUsageTable usage_;

        static const Enum<pointAction> pointActionNames_;

        //- Names of sets to use
        wordList names_;

        //- Option
        pointAction option_;


    // Private Member Functions

        //- Depending on point-to-cell option add to or delete from cellSet.
        void combine(topoSet& set, const bool add, const word& setName) const;


public:

    //- Runtime type information
    TypeName("pointToCell");

    // Constructors

        //- Construct from components
        pointToCell
        (
            const polyMesh& mesh,
            const word& setName,
            const pointAction option
        );

        //- Construct from dictionary
        pointToCell(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        pointToCell(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~pointToCell() = default;


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
