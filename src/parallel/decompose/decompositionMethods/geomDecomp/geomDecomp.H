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
    Foam::geomDecomp

Description
    Base for geometrical domain decomposition methods

    Base coefficients:
    \table
        Property  | Description                 | Required | Default
        n         | (nx ny nz)                  | yes
        delta     | delta for rotation matrix   | no | 0.001
    \endtable

SourceFiles
    geomDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef geomDecomp_H
#define geomDecomp_H

#include "decompositionMethod.H"
#include "Vector.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class geomDecomp Declaration
\*---------------------------------------------------------------------------*/

class geomDecomp
:
    public decompositionMethod
{
    // Private Member Functions

        //- Read input values and initialize the rotDelta_
        void readCoeffs();


protected:

    // Protected data

        //- Coefficients for all derived methods
        const dictionary& coeffsDict_;

        Vector<label> n_;

        //- Default = 0.001
        scalar delta_;

        tensor rotDelta_;

public:

    // Constructors

        //- Construct for derived type name and decomposition dictionary
        geomDecomp
        (
            const word& derivedType,
            const dictionary& decompDict,
            int select = selectionType::DEFAULT
        );

        //- Construct for derived type name, decomposition dictionary
        //- and region name
        geomDecomp
        (
            const word& derivedType,
            const dictionary& decompDict,
            const word& regionName,
            int select = selectionType::DEFAULT
        );


        //- Return for every coordinate the wanted processor number.
        virtual labelList decompose
        (
            const pointField& points,
            const scalarField& pointWeights
        ) = 0;

        //- Like decompose but with uniform weights on the points
        virtual labelList decompose(const pointField&) = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
