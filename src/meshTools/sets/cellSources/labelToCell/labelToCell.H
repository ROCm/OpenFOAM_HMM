/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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
    Foam::labelToCell

Description
    A topoSetCellSource to select cells based on explicitly given labels.

    \heading Dictionary parameters
    \table
        Property    | Description                       | Required  | Default
        value       | The cell labels                   | yes   |
    \endtable

SourceFiles
    labelToCell.C

\*---------------------------------------------------------------------------*/

#ifndef labelToCell_H
#define labelToCell_H

#include "topoSetCellSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class labelToCell Declaration
\*---------------------------------------------------------------------------*/

class labelToCell
:
    public topoSetCellSource
{

    // Private data

        //- Add usage string
        static addToUsageTable usage_;

        //- Cell labels read from dictionary
        labelList labels_;


public:

    //- Runtime type information
    TypeName("labelToCell");

    // Constructors

        //- Construct from components, copying labels
        labelToCell(const polyMesh& mesh, const labelList& labels);

        //- Construct from components, moving labels
        labelToCell(const polyMesh& mesh, labelList&& labels);

        //- Construct from dictionary
        labelToCell(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        labelToCell(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~labelToCell() = default;


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
