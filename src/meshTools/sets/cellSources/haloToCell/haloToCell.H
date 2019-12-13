/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::haloToCell

Description
    A topoSetCellSource to add/remove halo cells (cells attached to
    the outside of the current cellSet).
    When removing halo cells, this strips off any cells on the boundary
    of the set.

    \heading Dictionary parameters
    \table
        Property    | Description                           | Required | Default
        steps       | Number of grow/shrink steps to use    | no  | 1
    \endtable

SourceFiles
    haloToCell.C

\*---------------------------------------------------------------------------*/

#ifndef haloToCell_H
#define haloToCell_H

#include "topoSetCellSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class haloToCell Declaration
\*---------------------------------------------------------------------------*/

class haloToCell
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
    TypeName("haloToCell");


    // Constructors

        //- Construct from components
        explicit haloToCell(const polyMesh& mesh, const label steps=1);

        //- Construct from dictionary
        haloToCell(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        haloToCell(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~haloToCell() = default;


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
