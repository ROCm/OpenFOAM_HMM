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
    Foam::cellSelection

Description
    Cell selection methods in subsetMesh

SourceFiles
    cellSelection.C

\*---------------------------------------------------------------------------*/

#ifndef cellSelection_H
#define cellSelection_H

#include "dictionary.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;

/*---------------------------------------------------------------------------*\
                         Class cellSelection Declaration
\*---------------------------------------------------------------------------*/

class cellSelection
{
protected:

    // Protected data

        //- Name
        const word name_;

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Input dictionary
        const dictionary dict_;

public:

    //- Runtime type information
    TypeName("cellSelection");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            cellSelection,
            dictionary,
            (
                const word& name,
                const polyMesh& mesh,
                const dictionary& dict
            ),
            (name, mesh, dict)
        );


    // Constructors

        //- Construct from dictionary
        cellSelection
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


    // Selectors

        //- Return a reference to the selected cellSelection
        static autoPtr<cellSelection> New
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~cellSelection();


    // Member Functions

        //- Count global number of selected elements
        static label count(const boolList&);

        const word& name() const
        {
            return name_;
        }

        virtual void select(boolList&) const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
