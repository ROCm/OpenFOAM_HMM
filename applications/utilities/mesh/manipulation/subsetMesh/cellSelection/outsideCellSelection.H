/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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
    Foam::cellSelections::outsideCellSelection

Description
    Deselect cells not reachable from 'inside' points

SourceFiles
    outsideCellSelection.C

\*---------------------------------------------------------------------------*/

#ifndef outsideCellSelection_H
#define outsideCellSelection_H

#include "cellSelection.H"
#include "pointField.H"
#include "boolList.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class regionSplit;

namespace cellSelections
{

/*---------------------------------------------------------------------------*\
                     Class outsideCellSelection Declaration
\*---------------------------------------------------------------------------*/

class outsideCellSelection
:
    public cellSelection
{
    // Private data

        //- Locations to keep
        const pointField locationsInMesh_;

        //- Number of layers to erode
        const label nErode_;


   // Private Member Functions

        //- For debugging: generate volScalarField with 1.0 for all true
        tmp<volScalarField> generateField
        (
            const word& name,
            const boolList& lst
        ) const;

        //- Mark faces inbetween selected and unselected elements
        void markRegionFaces
        (
            const boolList& selectedCell,
            boolList& regionFace
        ) const;

        //- Determine for every disconnected region in the mesh whether
        //  it contains a locationInMesh
        boolList findRegions(const bool verbose, const regionSplit&) const;

        //- Unselect regions not containing a locationInMesh
        void unselectOutsideRegions(boolList& selectedCell) const;

        //- Unselect one layer of cells from selectedCell
        void shrinkRegions(boolList& selectedCell) const;

        //- Erode a given number of layers from selectedCell. Remove any
        //  region that gets disconnected that way.
        void erode(boolList& selectedCell) const;

public:

    //- Runtime type information
    TypeName("outside");


    // Constructors

        //- Construct from dictionary
        outsideCellSelection
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Clone
        autoPtr<cellSelection> clone() const
        {
            notImplemented("autoPtr<cellSelection> clone() const");
            return autoPtr<cellSelection>(NULL);
        }


    //- Destructor
    virtual ~outsideCellSelection();


    // Member Functions

        //- Apply this selector
        virtual void select(boolList&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace cellSelections
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
