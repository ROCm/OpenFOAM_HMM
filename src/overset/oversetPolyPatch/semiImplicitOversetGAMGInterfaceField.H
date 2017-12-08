/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::semiImplicitOversetGAMGInterfaceField

Description
    GAMG agglomerated processor interface field.

SourceFiles
    processorGAMGInterfaceField.C

\*---------------------------------------------------------------------------*/

#ifndef semiImplicitOversetGAMGInterfaceField_H
#define semiImplicitOversetGAMGInterfaceField_H

#include "oversetGAMGInterfaceField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class semiImplicitOversetGAMGInterfaceField Declaration
\*---------------------------------------------------------------------------*/

class semiImplicitOversetGAMGInterfaceField
:
    public oversetGAMGInterfaceField
{
    // Private Member Functions

        //- Disallow default bitwise copy construct
        semiImplicitOversetGAMGInterfaceField
        (
            const semiImplicitOversetGAMGInterfaceField&
        );

        //- Disallow default bitwise assignment
        void operator=(const semiImplicitOversetGAMGInterfaceField&);


public:

    //- Runtime type information
    TypeName("semiImplicitOverset");


    // Constructors

        //- Construct from GAMG interface and fine level interface field
        semiImplicitOversetGAMGInterfaceField
        (
            const GAMGInterface& GAMGCp,
            const lduInterfaceField& fineInterface
        );


    //- Destructor
    virtual ~semiImplicitOversetGAMGInterfaceField();


    // Member Functions

        //- Update result field based on interface functionality
        virtual void updateInterfaceMatrix
        (
            scalarField& result,
            const bool add,
            const scalarField& psiInternal,
            const scalarField& coeffs,
            const direction cmpt,
            const Pstream::commsTypes commsType
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
