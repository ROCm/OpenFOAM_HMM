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
    Foam::solidThermo

Description
    Fundamental solid thermodynamic properties

SourceFiles
    solidThermo.C

\*---------------------------------------------------------------------------*/

#ifndef solidThermo_H
#define solidThermo_H

#include "volFields.H"
#include "runTimeSelectionTables.H"
#include "fvMesh.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "basicSolidMixture.H"
#include "basicThermo.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class solidThermo Declaration
\*---------------------------------------------------------------------------*/

class solidThermo
:
    public basicThermo
{

protected:

    // Protected data

        //- Density field [kg/m^3]
        //  Named 'rhoThermo' to avoid (potential) conflict with solver density
        volScalarField rho_;


public:

    //- Runtime type information
    TypeName("solidThermo");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        solidThermo,
        fvMesh,
        (const fvMesh& mesh),
        (mesh)
    );

     // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        solidThermo,
        dictionary,
        (const fvMesh& mesh, const dictionary& dict),
        (mesh, dict)
    );


    // Constructors

        //- Construct from mesh
        solidThermo(const fvMesh&);

        //- Construct from mesh and dict
        solidThermo(const fvMesh&, const dictionary& dict);

        //- Return a pointer to a new solidThermo created from
        //  the solidThermophysicalProperties dictionary
        static autoPtr<solidThermo> New(const fvMesh&);

        //- Return a pointer to a new solidThermo created from
        //  local dictionary
        static autoPtr<solidThermo> New(const fvMesh&, const dictionary&);


    //- Destructor
    virtual ~solidThermo();


    // Member functions

        // Fields derived from thermodynamic state variables

            //- Density [kg/m^3]
            virtual tmp<volScalarField> rho() const;

            //- Return non-const access to the local density field [kg/m^3]
            virtual volScalarField& rho();

            //- Thermal conductivity [W/m/K]
            virtual tmp<volVectorField> Kappa() const = 0;

            //- Absorption coefficient [1/m]
            virtual tmp<volScalarField> kappaRad() const = 0;

            //- Scatter coefficient
            virtual tmp<volScalarField> sigmaS() const = 0;

            //- Emissivity coefficient [1/m]
            virtual tmp<volScalarField> emissivity() const = 0;


        // Per patch calculation


            //- Anisotropic thermal conductivity [W/m/K]
            virtual tmp<vectorField> Kappa
            (
                const label patchI
            ) const = 0;

             //- Absorption coefficient [1/m]
            virtual tmp<scalarField> kappaRad
            (
                const label patchI
            ) const = 0;

            //- Scatter coefficient
            virtual tmp<scalarField> sigmaS
            (
                const label patchI
            ) const = 0;

            //- Emissivity coefficient [1/m]
            virtual tmp<scalarField> emissivity
            (
                const label patchI
            ) const = 0;


    // I-O

        //- Read thermophysicalProperties dictionary
        virtual bool read();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
