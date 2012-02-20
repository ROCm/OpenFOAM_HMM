/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::boundaryFirstRenumber

Description
    Cuthill-McKee renumbering

SourceFiles
    boundaryFirstRenumber.C

\*---------------------------------------------------------------------------*/

#ifndef boundaryFirstRenumber_H
#define boundaryFirstRenumber_H

#include "renumberMethod.H"
#include "Switch.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class boundaryFirstRenumber Declaration
\*---------------------------------------------------------------------------*/

class boundaryFirstRenumber
:
    public renumberMethod
{
    // Private data

        const Switch reverse_;


    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        void operator=(const boundaryFirstRenumber&);
        boundaryFirstRenumber(const boundaryFirstRenumber&);


public:

    //- Runtime type information
    TypeName("boundaryFirst");


    // Constructors

        //- Construct given the renumber dictionary
        boundaryFirstRenumber(const dictionary& renumberDict);


    //- Destructor
    virtual ~boundaryFirstRenumber()
    {}


    // Member Functions

        //- Return the order in which cells need to be visited, i.e.
        //  from ordered back to original cell label.
        //  This is only defined for geometric renumberMethods.
        virtual labelList renumber(const pointField&) const
        {
            notImplemented
            (
                "boundaryFirstRenumber::renumber(const pointField&)"
            );
            return labelList(0);
        }

        //- Return the order in which cells need to be visited, i.e.
        //  from ordered back to original cell label.
        //  Use the mesh connectivity (if needed)
        virtual labelList renumber
        (
            const polyMesh& mesh,
            const pointField& cc
        ) const;

        //- Return the order in which cells need to be visited, i.e.
        //  from ordered back to original cell label.
        //  The connectivity is equal to mesh.cellCells() except
        //  - the connections are across coupled patches
        virtual labelList renumber
        (
            const labelListList& cellCells,
            const pointField& cc
        ) const
        {
            notImplemented
            (
                "boundaryFirstRenumber::renumber"
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
