/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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
    Foam::chemistrySolver

Description
    An abstract class for solving chemistry

SourceFiles
    chemistrySolver.C

\*---------------------------------------------------------------------------*/

#ifndef chemistrySolver_H
#define chemistrySolver_H

#include "ODEChemistryModel.H"
#include "IOdictionary.H"
#include "scalarField.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class chemistrySolver Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoType>
class chemistrySolver
{
protected:

    // Protected data

        //- Reference to the chemistry model
        ODEChemistryModel<ThermoType>& model_;


public:

        //- Runtime type information
        TypeName("chemistrySolver");


        // Declare runtime constructor selection table
        declareRunTimeSelectionTable
        (
            autoPtr,
            chemistrySolver,
            dictionary,
            (
                ODEChemistryModel<ThermoType>& model
            ),
            (model)
        );


    // Constructors

        //- Construct from components
        chemistrySolver(ODEChemistryModel<ThermoType>& model);


    //- Selector
    static autoPtr<chemistrySolver> New
    (
        ODEChemistryModel<ThermoType>& model
    );


    //- Destructor
    virtual ~chemistrySolver();


    // Member Functions

        //- Update the concentrations and return the chemical time
        virtual scalar solve
        (
            scalarField &c,
            const scalar T,
            const scalar p,
            const scalar t0,
            const scalar dt
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeChemistrySolver(Thermo)                                           \
                                                                              \
defineTemplateTypeNameAndDebug                                                \
(                                                                             \
    chemistrySolver<Thermo>,                                                  \
    0                                                                         \
);                                                                            \
                                                                              \
defineTemplateRunTimeSelectionTable(chemistrySolver<Thermo>, dictionary);


#define makeChemistrySolverType(SS, Thermo)                                   \
                                                                              \
    defineNamedTemplateTypeNameAndDebug(SS<Thermo>, 0);                       \
                                                                              \
    chemistrySolver<Thermo>::adddictionaryConstructorToTable<SS<Thermo> >     \
        add##SS##Thermo##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "chemistrySolver.C"
#   include "newChemistrySolver.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
