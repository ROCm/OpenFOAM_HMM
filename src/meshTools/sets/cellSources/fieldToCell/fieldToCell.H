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
    Foam::fieldToCell

Description
    A topoSetCellSource to select cells based on field values.

    \heading Dictionary parameters
    \table
        Property    | Description                       | Required  | Default
        field       | The (scalar, vector) field to use | yes   |
        min         | The min value for the subset      | yes   |
        max         | The max value for the subset      | yes   |
    \endtable

SourceFiles
    fieldToCell.C

\*---------------------------------------------------------------------------*/

#ifndef fieldToCell_H
#define fieldToCell_H

#include "topoSetCellSource.H"
#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fieldToCell Declaration
\*---------------------------------------------------------------------------*/

class fieldToCell
:
    public topoSetCellSource
{

    // Private data

        //- Add usage string
        static addToUsageTable usage_;


        //- Name of volScalarField, volVectorField
        word fieldName_;

        //- Min to subset field values with
        scalar min_;

        //- Max to subset field values with
        scalar max_;


    // Private Member Functions

        //- Depending on field values add to or delete from cellSet.
        void applyToSet
        (
            const topoSetSource::setAction action,
            const scalarField& field,
            topoSet& set
        ) const;

public:

    //- Runtime type information
    TypeName("fieldToCell");

    // Constructors

        //- Construct from components
        fieldToCell
        (
            const polyMesh& mesh,
            const word& fieldName,
            const scalar min,
            const scalar max
        );

        //- Construct from dictionary
        fieldToCell(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        fieldToCell(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~fieldToCell() = default;


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
