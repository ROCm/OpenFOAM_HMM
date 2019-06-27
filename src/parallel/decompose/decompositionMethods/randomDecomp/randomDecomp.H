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

        //- Random distribution on the 0-nCells interval
        labelList randomMap(const label nCells) const;

        //- No copy construct
        void operator=(const randomDecomp&) = delete;

        //- No copy assignment
        randomDecomp(const randomDecomp&) = delete;


public:

    //- Runtime type information
    TypeName("random");


    // Constructors

        //- Construct given the decomposition dictionary
        randomDecomp(const dictionary& decompDict);

        //- Construct for decomposition dictionary and region name
        randomDecomp
        (
            const dictionary& decompDict,
            const word& regionName
        );


    //- Destructor
    virtual ~randomDecomp() = default;


    // Member Functions

        //- Does not care about proc boundaries
        virtual bool parallelAware() const
        {
            return true;
        }

        //- Return for every coordinate the wanted processor number.
        virtual labelList decompose
        (
            const polyMesh& mesh,
            const pointField& cc,
            const scalarField& cWeights
        ) const;

        //- Return for every coordinate the wanted processor number.
        //  Explicitly provided connectivity - does not use mesh_.
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
