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
    Foam::directMappedVariableThicknessWallFvPatch

Description
    Take thickness field and number of layers and returns deltaCoeffs
    as 2.0/thickness/nLayers.
    To be used with 1D thermo baffle.

SourceFiles
    directMappedVariableThicknessWallFvPatch.C

\*---------------------------------------------------------------------------*/

#ifndef directMappedWallFvPatch_H
#define directMappedWallFvPatch_H

#include "wallFvPatch.H"
#include "directMappedVariableThicknessWallPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
        Class directMappedVariableThicknessWallFvPatch Declaration
\*---------------------------------------------------------------------------*/

class directMappedVariableThicknessWallFvPatch
:
    public wallFvPatch
{

protected:

    // Protected Member Functions

        //- Read neighbour cell distances from dictionary
        void makeDeltaCoeffs(scalarField& dc) const;


public:

    //- Runtime type information
    TypeName(directMappedVariableThicknessWallPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        directMappedVariableThicknessWallFvPatch
        (
            const polyPatch& patch,
            const fvBoundaryMesh& bm
        )
        :
            wallFvPatch(patch, bm)
        {}

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
