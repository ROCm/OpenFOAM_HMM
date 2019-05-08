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
    Foam::randomDecomp

Description
    Decomposition according to pseudo-random number generator.

SourceFiles
    randomDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef randomDecomp_H
#define randomDecomp_H

#include "decompositionMethod.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class randomDecomp Declaration
\*---------------------------------------------------------------------------*/

class randomDecomp
:
    public decompositionMethod
{
    // Private Member Functions

        //- Disallow default bitwise copy construct and assignment
        void operator=(const randomDecomp&);
        randomDecomp(const randomDecomp&);


public:

    //- Runtime type information
    TypeName("random");


    // Constructors

        //- Construct given the decomposition dictionary
        randomDecomp(const dictionary& decompositionDict);

    //- Destructor
    virtual ~randomDecomp()
    {}


    // Member Functions

        //- Manual decompose does not care about proc boundaries - is all
        //  up to the user.
        virtual bool parallelAware() const
        {
            return true;
        }

        //- Return for every coordinate the wanted processor number. Use the
        //  mesh connectivity (if needed)
        virtual labelList decompose
        (
            const polyMesh& mesh,
            const pointField& cc,
            const scalarField& cWeights
        ) const;

        //- Return for every coordinate the wanted processor number. Explicitly
        //  provided connectivity - does not use mesh_.
        //  The connectivity is equal to mesh.cellCells() except for
        //  - in parallel the cell numbers are global cell numbers (starting
        //    from 0 at processor0 and then incrementing all through the
        //    processors)
        //  - the connections are across coupled patches
        virtual labelList decompose
        (
            const labelListList& globalCellCells,
            const pointField& cc,
            const scalarField& cWeights
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
