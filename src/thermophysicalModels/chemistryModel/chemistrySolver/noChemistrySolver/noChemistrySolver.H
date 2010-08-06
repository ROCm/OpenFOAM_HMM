/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2010-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::noChemistry

Description
    Dummy chemistry solver for 'none' option

SourceFiles
    noChemistrySolver.H
    noChemistrySolver.C

\*---------------------------------------------------------------------------*/

#ifndef noChemistySolver_H
#define noChemistySolver_H

#include "chemistrySolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
template<class CompType, class ThermoType>
class noChemistrySolver;

/*---------------------------------------------------------------------------*\
                      Class noChemistrySolver Declaration
\*---------------------------------------------------------------------------*/

template<class CompType, class ThermoType>
class noChemistrySolver
:
    public chemistrySolver<CompType, ThermoType>
{

public:

    //- Runtime type information
    TypeName("none");


    // Constructors


        //- Construct from components
        noChemistrySolver
        (
            ODEChemistryModel<CompType, ThermoType>& model,
            const word& modelName
        );


    //- Destructor
    virtual ~noChemistrySolver();


    // Member Functions

        //- Update the concentrations and return the chemical time
        scalar solve
        (
            scalarField &c,
            const scalar T,
            const scalar p,
            const scalar t0,
            const scalar dt
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "noChemistrySolver.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
