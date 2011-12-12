/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::manualRenumber

Description
    Renumber given a cell-to-new cell association in a file

SourceFiles
    manualRenumber.C

\*---------------------------------------------------------------------------*/

#ifndef manualRenumber_H
#define manualRenumber_H

#include "renumberMethod.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class manualRenumber Declaration
\*---------------------------------------------------------------------------*/

class manualRenumber
:
    public renumberMethod
{
    // Private data

        const fileName dataFile_;


    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        void operator=(const manualRenumber&);
        manualRenumber(const manualRenumber&);


public:

    //- Runtime type information
    TypeName("manual");


    // Constructors

        //- Construct given the renumber dictionary
        manualRenumber(const dictionary& renumberDict);

    //- Destructor
    virtual ~manualRenumber()
    {}


    // Member Functions

        //- Return for every coordinate the wanted processor number.
        //  We need a polyMesh (to be able to load the file)
        virtual labelList renumber(const pointField&)
        {
            notImplemented("manualRenumber::renumber(const pointField&)");
            return labelList(0);
        }

        //- Return for every coordinate the wanted processor number. Use the
        //  mesh connectivity (if needed)
        virtual labelList renumber
        (
            const polyMesh& mesh,
            const pointField& cc
        );

        //- Return for every cell the new cell label.
        //  The connectivity is equal to mesh.cellCells() except
        //  - the connections are across coupled patches
        virtual labelList renumber
        (
            const labelListList& cellCells,
            const pointField& cc
        )
        {
            notImplemented
            (
                "manualRenumber::renumber"
                "(const labelListList&, const pointField&)"
            );
            return labelList(0);
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
