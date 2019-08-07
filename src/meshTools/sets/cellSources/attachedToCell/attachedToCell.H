/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
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
    Foam::attachedToCell

Description
    A topoSetCellSource to add/remove cells attached to the current cellSet

    \heading Dictionary parameters
    \table
        Property    | Description                           | Required | Default
        steps       | Number of grow/shrink steps to use    | no  | 1
    \endtable

SourceFiles
    attachedToCell.C

\*---------------------------------------------------------------------------*/

#ifndef attachedToCell_H
#define attachedToCell_H

#include "topoSetCellSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class attachedToCell Declaration
\*---------------------------------------------------------------------------*/

class attachedToCell
:
    public topoSetCellSource
{

    // Private data

        //- Add usage string
        static addToUsageTable usage_;

        //- Number of grow/shrink steps
        label steps_;


    // Private Member Functions

        void combine(topoSet& set, const bool add) const;


public:

    //- Runtime type information
    TypeName("attachedToCell");


    // Constructors

        //- Construct from components
        explicit attachedToCell(const polyMesh& mesh, const label steps=1);

        //- Construct from dictionary
        attachedToCell(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        attachedToCell(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~attachedToCell() = default;


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
