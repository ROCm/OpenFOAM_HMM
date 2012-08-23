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
    Foam::basicThermo

Description
    Basic thermodynamic properties

SourceFiles
    basicThermo.C
    basicThermoNew.C

\*---------------------------------------------------------------------------*/

#ifndef basicThermo_H
#define basicThermo_H

#include "veryBasicThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class basicThermo Declaration
\*---------------------------------------------------------------------------*/

class basicThermo
:
    public veryBasicThermo
{

public:

    //- Runtime type information
    TypeName("basicThermo");


    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        basicThermo,
        fvMesh,
        (const fvMesh& mesh),
        (mesh)
    );

    // Constructors

        //- Construct from mesh
        basicThermo(const fvMesh&);

        //- Construct from mesh
        basicThermo(const fvMesh&, const dictionary&);


    //- Selector
    static autoPtr<basicThermo> New(const fvMesh&);

    //- Destructor
    virtual ~basicThermo();


    // Member functions

        // Access to thermodynamic state variables

            //- Compressibility [s^2/m^2]
            virtual const volScalarField& psi() const = 0;


        // Access to transport state variables

            //- Dynamic viscosity of mixture [kg/m/s]
            virtual const volScalarField& mu() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
