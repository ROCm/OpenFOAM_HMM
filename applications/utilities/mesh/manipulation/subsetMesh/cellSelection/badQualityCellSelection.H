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
    Foam::cellSelections::badQualityCellSelection

Description
    Deselect bad quality cells

SourceFiles
    badQualityCellSelection.C

\*---------------------------------------------------------------------------*/

#ifndef badQualityCellSelection_H
#define badQualityCellSelection_H

#include "cellSelection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace cellSelections
{

/*---------------------------------------------------------------------------*\
                     Class badQualityCellSelection Declaration
\*---------------------------------------------------------------------------*/

class badQualityCellSelection
:
    public cellSelection
{

public:

    //- Runtime type information
    TypeName("badQuality");


    // Constructors

        //- Construct from dictionary
        badQualityCellSelection
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
    virtual ~badQualityCellSelection();


    // Member Functions

        virtual void select(boolList&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace cellSelections
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
